package SysBio;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

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
	public static void main(String ar[]) throws IOException, IloException
	{
		/*---------Creating Genes-----------------*/
		HashMap<String, Gene> Genes = new HashMap<String, Gene>();
		File nodes = new File("H:\\Beethika\\2nd sem\\Systems Biology\\Project\\newNodes.txt");
		BufferedReader br = new BufferedReader(new FileReader(nodes));
		String s = br.readLine();
		while(br.ready())
		{
			s = br.readLine();
			String[] l = s.split("\\t");
			Gene g = new Gene();
			g.GeneId = l[0].split("\\|")[1];
			//System.out.println(g.GeneId);
			g.pval = Double.parseDouble(l[1]);
			g.Regulation = Integer.parseInt(l[2]);
			g.degree = 0;
			g.costScore = 0;
			g.profitScore = 0;
			Genes.put(g.GeneId, g);
		}
		br.close();
		
		/*---------mapping Proteins to Genes-----------------*/
		HashMap<String, Gene> proToGene = new HashMap<String, Gene>();
		File gene2ensembl = new File("H:\\Beethika\\2nd sem\\Systems Biology\\Project\\gene2ensembl_dummy.txt");
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
		File inFile = new File("H:\\Beethika\\2nd sem\\Systems Biology\\Project\\data\\pro_action_dummy.txt");      //Reading from file
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
		// --- verify for the threshold and constant???????????
    	double threshold = 0.05; 
    	double cons = Math.log(0.05);
    	double maxWn = 0;
    	double maxCn = 0;
    	/*for(Gene g :G)
    		System.out.println(g);*/
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
    		if(g.pval < threshold) // significant genes
    			g.profitScore = -1*Math.log(g.pval);
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
    	
    	/* computing the edge score */
    	double maxSe = 0;
    	for(Edge e : E)
    	{
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
    			
    			e.edgeScore = new ArrayList<Double>();
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
    	
    	/*-------- Writing  Gene Gene Interaction Network to File-------*/
    	FileWriter writer1 = new FileWriter("H:\\Beethika\\2nd sem\\Systems Biology\\Project\\data\\CompleteGGI.txt", false);
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
    	writer2.close();
    	
    	ASN(E,G);
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
	public static void ASN(Collection<Edge> E,Collection<Gene> G) throws IloException
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
			System.out.println(g.GeneId);
			c++;
		}
		
		//double Se[] = new;
		c=0;
		ArrayList<Double> Se=new ArrayList<Double>();
		ArrayList<Integer> xe_n1=new ArrayList<Integer>();
		ArrayList<Integer> xe_n2=new ArrayList<Integer>();
		
		for(Edge e : E)
		{
			for(int i=0;i<e.confidenceScore.size();i++)
			{
				Se.add(e.confidenceScore.get(i));
				if(e.Direction.get(i)==-1)
				{
					xe_n1.add(GeneMap.get(e.G2));
					xe_n2.add(GeneMap.get(e.G1));
					System.out.println(e.G2+" "+e.G1);
				}
				else 
				{
					xe_n1.add(GeneMap.get(e.G1));
					xe_n2.add(GeneMap.get(e.G2));
					System.out.println(e.G1+" "+e.G2);
				}
				c++;
			}
		}
		
		IloCplex model = new IloCplex();
		int n= G.size();
		int e =E.size();
		//variables
		IloNumVar xe[] = model.boolVarArray(e) ;
		IloNumVar xn[] = model.boolVarArray(n) ;
		// objective 
		IloLinearNumExpr obj = model.linearNumExpr();
		for ( int i=0;i<n;i++)
		{
			obj.addTerm(W_bC[i],xn[i]);
		}
		for (int i=0; i<e; i++)
		{
			obj.addTerm(Se.get(i), xe[i]);
		}
		model.addMaximize(obj);
		//constraints
		for (int i =0; i<n; i++)
		{
			//int tot_xe = 0;
			IloLinearNumExpr expr = model.linearNumExpr();
			for (int j=0; j<e; j++)
			{
				if( i == xe_n1.get(j))
				{
				expr.addTerm(1, xe[j]);
				}
			}
			model.addLe(xn[i], expr);
		}
		for (int i=0; i<e; i++)
		{
			IloLinearNumExpr expr = model.linearNumExpr();
			for (int j=0; j<n; j++)
			{
				if (xe_n1.get(i) == j)
				{
				expr.addTerm(1, xn[j]);
				model.addLe(xe[i], expr);
				break;
				}
			}
			
		}
		for (int i=0; i<e; i++)
		{
			IloLinearNumExpr expr = model.linearNumExpr();
			for (int j=0; j<n; j++)
			{
				if (xe_n2.get(i) == j)
				{
					expr.addTerm(1, xn[j]);
					model.addLe(xe[i], expr);
					break;
				}
			}
			
		}
		model.solve();
		System.out.println(" no of constraints "+model.getNrows());
		model.output().println("Solution status = " + model.getStatus());
		model.output().println("Solution value = " + model.getObjValue());
		for(int i=0;i<xe.length;i++)
			System.out.println("edges"+model.getValue(xe[i]));
		for(int i=0;i<xn.length;i++)
			System.out.println("nodes"+model.getValue(xn[i]));
		model.end();
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
	double pval;
	int Regulation;
	int degree;
	double profitScore;
	double costScore;
}