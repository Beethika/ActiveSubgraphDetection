package Gurobi;



import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

public class PPI
{
	public static void main(String ar[]) throws IOException, GRBException
	{
		/*---------Creating Genes-----------------*/
		HashMap<String, Gene> Genes = new HashMap<String, Gene>();
		File nodes = new File("H:\\Beethika\\2nd sem\\Systems Biology\\Project\\nodes.txt");
		BufferedReader br = new BufferedReader(new FileReader(nodes));
		String s = br.readLine();
		while(br.ready())
		{
			s = br.readLine();
			String[] l = s.split("\\t");
			Gene g = new Gene();
			g.GeneId = l[0].split("\\|")[1];
			g.GeneName = l[0].split("\\|")[0];
			//System.out.println(g.GeneId);
			g.pval = Double.parseDouble(l[1]);
			g.Regulation = Integer.parseInt(l[2]);
			g.degree = 0;
			g.costScore = 0;
			g.profitScore = 0;
			//if(g.Regulation!=0)
				Genes.put(g.GeneId, g);
		}
		br.close();
		
		/*---------mapping Proteins to Genes-----------------*/
		HashMap<String, Gene> proToGene = new HashMap<String, Gene>();
		File gene2ensembl = new File("H:\\Beethika\\2nd sem\\Systems Biology\\Project\\gene2ensembl");
		br = new BufferedReader(new FileReader(gene2ensembl));
		s = br.readLine();
		while(br.ready())
		{
			s = br.readLine().trim();
    		String[] l = s.split("\\t");
    		String geneId = l[1];
    		String protiene = l[6];
    		if(protiene.length()>1&&(Genes.containsKey(geneId) && !proToGene.containsKey(protiene)))
    		{
    			//System.out.println(Genes.containsKey(geneId));
    			//System.out.print(geneId+" "+protiene+" "+protiene.length()+"\n");
    			proToGene.put(protiene, Genes.get(geneId));
    		}
		}
		br.close();
		System.out.println("pro  to gene "+proToGene.size()+" Genes "+Genes.size());
		/*------------------ Reading PPI and converting to GGI-----------------*/ 
		//9606.protein.actions.v10.txt
		File inFile = new File("H:\\Beethika\\2nd sem\\Systems Biology\\Project\\data\\9606.protein.actions.v10.txt");      //Reading from file
    	br = new BufferedReader(new FileReader(inFile));
    	
    	HashMap<String,Integer> action= new HashMap<String, Integer>();
		action.put("binding", 1);
		action.put("activation", 2);
		action.put("inhibition",3);
		action.put("catalysis", 4);
		action.put("ptmod", 5);
		
		//action.put("reaction", 6); //////////???????
		//action.put("expression", 7); //////////???????????
		HashMap<String,Edge> protieneMap = new HashMap<String, Edge>();
    	s = br.readLine();
    	double maxCe = 0;
    	while (br.ready())
    	{
    		s = br.readLine();
    		String[] l = s.split("\\t");
    		if((l[2].equals("activation")&&l[4].equals("0"))||(l[2].equals("inhibition")&&l[4].equals("0"))||l[2].equals("expression")||l[2].equals("reaction"))
    			continue;
    		else
    			maxCe = addEdge(l,protieneMap,action,proToGene,maxCe);	
    	}
    	br.close();
    	System.out.println("protiene map "+ protieneMap.size());
    	/* Standardising confidence score of edge */
       	Collection<Edge> E =  protieneMap.values();
       	for(Edge e : E)
       	{
       		for(int i =0;i<e.confidenceScore.size();i++)
       		{
       			Double ce = e.confidenceScore.remove(i);
    			e.confidenceScore.add(i,ce / maxCe);
       		}
       	}
       	
    	/* Genes in Protein to Gene Mapping are the only Genes present in the Network
    	 * So construct the Node Score For only these Genes
    	 */
    	Collection<Gene> G =  proToGene.values();
    	//System.out.println(G.isEmpty());
		
    	/*for(Gene g :G)
    		System.out.println(g);*/
    	
    	
    	/* computing the edge score */
    	double maxSe = 0;
    	for(Edge e : E)
    	{
    		e.edgeScore = new ArrayList<Double>();
    		for(int i = 0; i < e.Action.size() ; i++)
    		{
    			Gene g1 = Genes.get(e.G1);
    			Gene g2 = Genes.get(e.G2);
    			int t = e.Action.get(i);
    			double c = e.confidenceScore.get(i);
    			double w = 0;
    			if(t == 2)
    			{
    				if((g1.Regulation == 1 && g2.Regulation == 1) || (g1.Regulation == -1 && g2.Regulation == -1))
    					w = 2;
    				else if (g1.Regulation == 0 || g2.Regulation == 0)
    					w = -1;
    				else if((g1.Regulation == 1 && g2.Regulation == -1) || (g1.Regulation == -1 && g2.Regulation == 1))
    					w = -2;
    			}
    			else if(t == 3)
    			{
    				if((g1.Regulation == 1 && g2.Regulation == 1) || (g1.Regulation == -1 && g2.Regulation == -1))
    					w = -2;
    				else if (g1.Regulation == 0 || g2.Regulation == 0)
    					w = -1;
    				else if((g1.Regulation == 1 && g2.Regulation == -1) || (g1.Regulation == -1 && g2.Regulation == 1))
    					w = 2;
    			}
    			else // for catalytic, binding and post translational modification
    			{
    				w = -1;
    			}
    			
    			
    			double se = w * c;
    			e.edgeScore.add(se);
    			if(maxSe < Math.abs(se))
    				maxSe = Math.abs(se);	
    		}
    	}
    	
    	/* Standardising Edge Score */
    	for(Edge e : E)
    	{
    		for(int i = 0; i < e.edgeScore.size() ; i++)
    		{
    			Double se = e.edgeScore.remove(i);
    			e.edgeScore.add(i,se / maxSe);
    		}
    	}
    	
    	double threshold[] = {0.05,0.005};//,0.01,0.001}; 
    	for(int i=0 ; i<threshold.length ;i++)
    	{
    		// --- verify for the threshold and constant???????????
    		//threshold /= 10;
        	//double tau = Math.pow(2,i);
    		
        	double cons = 1*Math.log(threshold[i]);
        	double maxWn = 0;
        	double maxCn = 0;
        	ArrayList<Gene> FinalGene = new ArrayList<Gene>();
    	
    	for(Gene g : G)
    	{
    		/* computing Profit */
    		// distinguish between the significant genes and not significant genes 
    		//System.out.println(g);
    		if(g.degree == 0)
    		{
    			//G.remove(g);
    			continue;
    		}
    		if(g.pval < threshold[i]) // significant genes
    			g.profitScore = -1*(Math.log(g.pval));//-Math.log(tau));
    		else
    			g.profitScore = cons;
    		if(maxWn < Math.abs(g.profitScore))
    			maxWn = Math.abs(g.profitScore);
    		
    		/* penalise for the highly connected node */
    		
    		g.costScore = Math.log(g.degree);
    		//System.out.println(g.degree+" "+g.costScore);
    		if(maxCn < Math.abs(g.costScore))
    			maxCn = Math.abs(g.costScore);
    		FinalGene.add(g);
		}
    	G = FinalGene;
    	System.out.println("gene size "+ G.size());
    	/* standardising the Node Score */
    	for(Gene g : G)
    	{
    		g.profitScore = g.profitScore / maxWn ;
    		g.costScore = g.costScore / maxCn ;
    	}
    	
    	/*-------- Writing  Gene Gene Interaction Network to File-------*/
    	/**FileWriter writer1 = new FileWriter("H:\\Beethika\\2nd sem\\Systems Biology\\Project\\data\\CompleteGGI.txt", false);
    	FileWriter writer2 = new FileWriter("H:\\Beethika\\2nd sem\\Systems Biology\\Project\\data\\OnlyGGI.txt", false);
     	//System.out.println(E.size());
    	for(Edge e: E)
    	{
    		writer1.write("Genes: \t"+e.G1+"\t"+e.G2+"\r\n");
    		writer1.write("Profit Score: \t"+Genes.get(e.G1).profitScore+"\t"+Genes.get(e.G2).profitScore+"\r\n");
    		writer1.write("Cost Score: \t"+Genes.get(e.G1).costScore+"\t"+Genes.get(e.G2).costScore+"\r\n");
    		writer2.write(e.G1+"\t"+e.G2+"\r\n");
    		writer1.write("Action: \t");
    		for(int t: e.Action)
 	   			writer1.write(t+"\t"); 
    		writer1.write("\r\n");
    		writer1.write("Confidence Score: \t");
    		for(double ce: e.confidenceScore)
 	   			writer1.write(ce+"\t");
    		writer1.write("\r\n");
    		writer1.write("Direction: \t");
    		for(double d: e.Direction)
 	   			writer1.write(d+"\t"); 
    		writer1.write("\r\n");
    		writer1.write("Edge Score: \t");
    		for(double Se: e.edgeScore)
 	   			writer1.write(Se+"\t"); 
    		writer1.write("\r\n");
    		writer1.write("\r\n");
    	}
    	writer1.close();
    	writer2.close();*/
    	
    	ASN(E,G,threshold[i]);
    	}
    	/* computing sub-network score 
    	 * done using MILP
    	 * Find Subnetwork based on maximising the sub-network score
    	 */
    	
	}
	static double addEdge(String l[],HashMap<String,Edge> protieneMap,HashMap<String,Integer> action,HashMap<String, Gene> proToGene,double maxCe)
	{
		String tag;
		int check = 0;
		String p1 = l[0].split("\\.")[1];
		String p2 = l[1].split("\\.")[1];
		//System.out.println(p1+" "+p2);
		if(l[0].compareTo(l[1])<0)
			tag = p1+"#"+p2;
		else
		{
			tag = p2+"#"+p1;
			check = 1;
		}
		if(proToGene.containsKey(p1)&&proToGene.containsKey(p2))// removing the proteins whose mapping not found
		{

			Edge e;
			//System.out.println("hello");
			if(protieneMap.containsKey(tag))
			{
				e = protieneMap.get(tag);
				if(!e.Action.contains(action.get(l[2])))
				{
					e.Action.add(action.get(l[2]));
					e.confidenceScore.add(Double.parseDouble(l[5]));
				}
			}
			else
			{
				e = new Edge();
				if(p1.compareTo(p2)<0)
				{
					e.G1 = proToGene.get(p1).GeneId;
					e.G2 = proToGene.get(p2).GeneId;
				}
				else
				{	
					e.G2 = proToGene.get(p1).GeneId;
					e.G1 = proToGene.get(p2).GeneId;
				}
				e.Action = new ArrayList<Integer>();
				e.Action.add(action.get(l[2]));
				e.confidenceScore = new ArrayList<Double>();
				e.confidenceScore.add(Double.parseDouble(l[5]));
				e.Direction = new ArrayList<Integer>();
				protieneMap.put(tag,e);
			}
			
			/* Adding Direction to the edge for different actions 
			 * 	Updating the degree of the nodes
			 */
			if(action.get(l[2]) == 2 || action.get(l[2]) == 3)
			{
				if(check == 0)
				{
					e.Direction.add(1);
				}
				else
				{
					e.Direction.add(-1);
				}
			}
			else
				e.Direction.add(0);
			Gene g1 = proToGene.get(p1);
			Gene g2 = proToGene.get(p2);
			g1.degree = g1.degree+1;
			g2.degree = g2.degree+1;
			if(maxCe < Math.abs(Double.parseDouble(l[5])))
				maxCe = Math.abs(Double.parseDouble(l[5]));
		}
		return maxCe;
	}
	public static void ASN(Collection<Edge> E,Collection<Gene> G,double tau) throws GRBException, IOException 
	{
		// data
		double[] W_bC = new double[G.size()];
		HashMap<String, Integer> GeneMap = new HashMap<String, Integer>();
		int c = 0;
		for(Gene g: G)
		{
			int bn = 0 ;
			if(g.profitScore<0)
			{
				bn = 1; 
			}
			W_bC[c] = g.profitScore - bn * g.costScore;
			GeneMap.put(g.GeneId,c);
			//System.out.println(g.GeneId);
			c++;
		}
		
		//double Se[] = new;
		c=0;
		ArrayList<Double> Se=new ArrayList<Double>();
		ArrayList<Integer> xe_n1=new ArrayList<Integer>();
		ArrayList<Integer> xe_n2=new ArrayList<Integer>();
		ArrayList<String> action = new ArrayList<String>();
		for(Edge e : E)
		{
			for(int i=0;i<e.edgeScore.size();i++)
			{
				Se.add(e.edgeScore.get(i));
				if(e.Direction.get(i)==-1)
				{
					xe_n1.add(GeneMap.get(e.G2));
					xe_n2.add(GeneMap.get(e.G1));
					//System.out.println(e.G2+" "+e.G1);
				}
				else 
				{
					xe_n1.add(GeneMap.get(e.G1));
					xe_n2.add(GeneMap.get(e.G2));
					//System.out.println(e.G1+" "+e.G2);
				}
				switch(e.Action.get(i))
				{
					case 1: action.add("binding"); break;
					case 2:	action.add("activation"); break;
					case 3: action.add("inhibition"); break;
					case 4: action.add("catalysis"); break;
					case 5: action.add("ptmod"); break;
				}
				c++;
			}
		}
		//System.load("H:\\gurobi651\\win64\\bin\\GurobiJni65.dll");
		GRBEnv env = new GRBEnv();
		GRBModel model = new GRBModel(env);
		
		int n= G.size();
		int e =xe_n1.size();
		//variables
		GRBVar xe[] = new GRBVar[e];
		for(int i=0;i<e;i++)
		 xe[i] = model.addVar(0,1,1, GRB.BINARY, "edges"+i) ;
		GRBVar xn[] = new GRBVar[n];
		for(int i=0;i<n;i++)
		 xn[i] = model.addVar(0,1,1, GRB.BINARY, "nodes"+i) ;
		model.update();
		// objective 
		GRBLinExpr obj = new GRBLinExpr();
		for ( int i=0;i<n;i++)
		{
			obj.addTerm(W_bC[i],xn[i]);
		}
		for (int i=0; i<e; i++)
		{
			obj.addTerm(Se.get(i), xe[i]);
		}
		model.setObjective(obj, GRB.MAXIMIZE);
		ArrayList<Gene> G1 = new ArrayList<Gene>(G);
		//constraints
		for (int i =0; i<n; i++)
		{
			//int tot_xe = 0;
			GRBLinExpr expr = new GRBLinExpr();
			for (int j=0; j<e; j++)
			{
				if( i == xe_n1.get(j))
				{
				expr.addTerm(1, xe[j]);
				}
			}
			model.addConstr(xn[i], GRB.LESS_EQUAL, expr, "NodeConstr"+i);
		}
		for (int i=0; i<e; i++)
		{
			GRBLinExpr expr = new GRBLinExpr();
			for (int j=0; j<n; j++)
			{
				if (xe_n1.get(i)==j)
				{
				expr.addTerm(1, xn[j]);
				model.addConstr(xe[i], GRB.LESS_EQUAL, expr, "EdgeN1Constr"+i);
				break;
				}
			}
			
		}
		for (int i=0; i<e; i++)
		{
			GRBLinExpr expr = new GRBLinExpr();
			for (int j=0; j<n; j++)
			{
				if (xe_n2.get(i)==j)
				{
					expr.addTerm(1, xn[j]);
					model.addConstr(xe[i], GRB.LESS_EQUAL, expr, "EdgeN2Constr"+i);
					break;
				}
			}
			
		}
		model.optimize();
		//System.out.println(" no of constraints "+model.getNrows());
		//model.output().println("Solution status = " + model.getStatus());
		//model.output().println("Solution value = " + model.getObjValue());
		String path = "H:\\Beethika\\2nd sem\\Systems Biology\\Project\\EdgesOutput_cons_new"+tau+".txt";
		FileWriter writer1 = new FileWriter(path, false);
		for(int i=0;i<xe.length;i++)
		{
			double o = xe[i].get(GRB.DoubleAttr.X);
			if(o == 1)
			{
				Gene g1 = G1.get(xe_n1.get(i));
				Gene g2 = G1.get(xe_n2.get(i));
				writer1.write(g1.GeneId+"."+g1.GeneName+"\t"+action.get(i)+"\t"+g2.GeneId+"."+g2.GeneName);
				writer1.write("\r\n");
			}
		}
		writer1.close();
		path = "H:\\Beethika\\2nd sem\\Systems Biology\\Project\\NodesOutput_cons_new"+tau+".txt";
		FileWriter writer2 = new FileWriter(path, false);
		c=0;
		ArrayList<Gene> ActiveGene = new ArrayList<Gene>();
		for(int i=0;i<xn.length;i++)
		{	
			double o = xn[i].get(GRB.DoubleAttr.X);
			if(o == 1)
			{
				ActiveGene.add(G1.get(i));
				writer2.write(G1.get(i).GeneId+"."+G1.get(i).GeneName+"\t"+G1.get(i).Regulation+"\t"+W_bC[c++]);
				writer2.write("\r\n");
			}
		}
		writer2.close();
		model.terminate();
		findingDrugTargets(ActiveGene);
	}
	
	//finding out Drug targets
	static void findingDrugTargets(ArrayList<Gene> ActiveGene) throws IOException
	{
		BufferedReader br = new BufferedReader(new FileReader(new File("H:\\Beethika\\2nd sem\\Systems Biology\\Project\\data\\Drugs.txt")));
		HashMap<String, String> Drug_Gene = new HashMap<String, String>();
		br.readLine();
		while(br.ready())
		{
			String l[] = br.readLine().split("\\t");
			String Value = l[0];
			String k = l[3];
			if(k.charAt(0)=='"')
				k=k.substring(1,k.length()-1);
				
			String s[] = k.split("\\,");
			for(int i=0;i<s.length;i++)
			{
				String key = s[i].trim();
				//System.out.println(key);
				Drug_Gene.put(key, Value);
			}
			
		}
		br.close();
		for(int i=0;i<ActiveGene.size();i++)
		{
			Gene g = ActiveGene.get((i));
			if(Drug_Gene.containsKey(g.GeneName))
			{
				System.out.println("For Active Gene "+g.GeneId+"."+g.GeneName+" Drug is "+Drug_Gene.get(g.GeneName));
			}
		}
		
	}
	
}

class Edge
{
	String G1;
	String G2;
	ArrayList<Integer> Action;
	ArrayList<Integer> Direction; // 1 for ->,-1 for <- and 0 for - 
	ArrayList<Double> confidenceScore;
	ArrayList<Double> edgeScore;
}
class Gene
{
	String GeneId;
	String GeneName;
	double pval;
	int Regulation;
	int degree;
	double profitScore;
	double costScore;
}