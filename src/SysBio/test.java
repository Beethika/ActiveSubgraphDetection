package SysBio;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class test {

	public static void main(String[] args) throws IOException 
	{
		File nodes = new File("H:\\Beethika\\2nd sem\\Systems Biology\\Project\\test.txt");
		BufferedReader br = new BufferedReader(new FileReader(nodes));
		String s = br.readLine();
		while(br.ready())
		{
			s = br.readLine();
			String[] l = s.split("\\t");
			System.out.println(l[0]+" "+l[1]);
			String k[] = l[0].split("\\.");
			for(String r : k)
				System.out.print(r+" ");
			System.out.println();
			
			//Genes.put(g.GeneId, g);
		}
		br.close();

	}

}
