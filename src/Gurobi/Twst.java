package Gurobi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class Twst
{
	public static void main(String ar[]) throws IOException
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
				System.out.println(key);
				Drug_Gene.put(key, Value);
			}
			
		}
		br.close();
		/*for(int i=0;i<ActiveGene.size();i++)
		{
			Gene g = ActiveGene.get((i));
			if(Drug_Gene.containsKey(g.GeneName))
			{
				System.out.println("For Active Gene "+g.GeneId+"."+g.GeneName+" Drug is "+Drug_Gene.get(g.GeneName));
			}
		}*/
		
	}

}
