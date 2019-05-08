import java.io.*;
import java.util.StringTokenizer;
import java.util.Hashtable;
import java.util.zip.GZIPInputStream; 
import java.math.BigDecimal; 
import java.math.*;


public class CountCurrent {


	static final int onset = 18;
	static final int nucleus = 17; // 17 == current; 
	static final int coda = 7;

	public static void main (String[] args) throws IOException {


		String fileName = args[0];
		//	BufferedReader d = new BufferedReader(new InputStreamReader(new FileInputStream(new File (fileName))));
		//                OutputStreamWriter out = new OutputStreamWriter (new FileOutputStream(fileName+".1"), "UTF-8");
		BufferedReader d;
		d = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileName))));
		String str = new String();
		str = d.readLine();

		int[][] ar = new int[onset][nucleus]; 
		for (int i=0; i<onset; i++) {
			for (int j=0; j<nucleus; j++) { ar[i][j] = 0; 
			}
		}
		int c = -1; int v=-1; 

		/* counting cv  */
		while (str != null) {
			str = str.trim();

			for (int i=0; i<str.length(); i++) {
				String tok = str.substring(i, i+1);  
				//System.out.println(tok); 
				if (i==0) { c = cNum(tok); }
				else if (i==1) { v = vNum(tok); }
			}
			int count = ar[c][v]; 
			count++; 
			ar[c][v] = count; 
			c = -1; v = -1; 

			str = d.readLine();
		}

		/* counting CV */ 
		Hashtable<String, Integer> hashCV = new Hashtable<String, Integer>();

		int sum = 0;
		int sumsum = 0; 

		for (int i=0; i<onset; i++) { // # of C;
			sum = 0; 
			for (int j=0; j<nucleus; j++) { // # of V; 
				sum +=  ar[i][j]; 
				sumsum += ar[i][j]; 
			}
			String iString = cString(i); 
			hashCV.put(iString, sum); 
			//System.out.println(iString + ":" + sum); 
		}
		sum = 0; 
		for (int j=0; j<nucleus; j++) { // # of V; 
			sum = 0; 
			for (int i=0; i<onset; i++) { // # of C;
				sum +=  ar[i][j]; 
			}
			String jString = vString(j); 
			hashCV.put(jString, sum); 
			//System.out.println(jString + ":" + sum); 

		}

		//System.out.println(sumsum); 


		// printing  c v numbers; 
		// 0 0 = # of ㄱ ㅏ 
		double rphi6CVsum = 0; 
		double ficherCVsum = 0; 
		double chiCVsum = 0; 
		int rphi6CVnum = 0; 
		double[] numArrayCV = new double[onset*nucleus]; 
		double[] ficherArrayCV = new double[onset*nucleus];
		double[] chiArrayCV = new double[onset*nucleus];

		for (int i=0; i<onset; i++) { // # of C;
			for (int j=0; j<nucleus; j++) { // # of V; 
				String iString = cString(i);
				String jString = vString(j);
				// printing all possible c v & their a b c d; 
				// a = # cv 
				// b = # cV - cv where v is the current v; 
				// c = # Cv - cv where c is the current c; 
				// d = # CV - cv where c/v are the current c/v; 

				int anum = ar[i][j]; 
				int bnum = hashCV.get(iString) - anum; 
				int cnum = hashCV.get(jString) - anum; 
				int dnum = sumsum - (anum+bnum+cnum); 

				double tpp = (double) anum / (anum + cnum);
				double deltapp = ((double) anum / (anum + cnum)) - ( (double) bnum / (bnum + dnum)); 
				double deltap = ((double) anum / (anum + bnum)) - ( (double) cnum / (cnum + dnum));

				double rphi5 = (double) ((double) (anum*dnum) - (bnum*cnum)) / Math.sqrt((anum + bnum) * (cnum + dnum) * (anum + cnum) * (bnum + dnum)); 
				double rphi6 = (double) Math.sqrt( deltap * deltapp);  


				FisherExact fe = new FisherExact(anum+bnum+cnum+dnum);
				double ficherCVdouble = fe.getCumlativeP(anum,bnum,cnum,dnum);



				if (rphi6<0) rphi6 = rphi6 * -1; 

				System.out.println( iString + "" + jString + "\t" + + rphi6 + "\t" + ficherCVdouble); 
				if (!Double.isNaN(rphi6)) numArrayCV[rphi6CVnum] =  rphi6; 
				else numArrayCV[rphi6CVnum] =  0; 
				if (!Double.isNaN(ficherCVdouble)) ficherArrayCV[rphi6CVnum] =  ficherCVdouble;
				else ficherArrayCV[rphi6CVnum] =  0;

				rphi6CVnum++; 
				if (!Double.isNaN(rphi6)) rphi6CVsum += rphi6;
				if (!Double.isNaN(ficherCVdouble)) ficherCVsum += ficherCVdouble;

			}
		}
		System.out.println(); 

		/* **************************************************************************************************************** */
		/* **************************************************************************************************************** */
		/* **************************************************************************************************************** */
		/* **************************************************************************************************************** */
		/* **************************************************************************************************************** */

		d = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileName))));
		str = d.readLine();

		ar = new int[nucleus][coda]; 
		for (int i=0; i<nucleus; i++) {
			for (int j=0; j<coda; j++) { ar[i][j] = 0; 
			}
		}
		c = -1; v=-1; 


		/* counting vc  */
		while (str != null) {
			str = str.trim();

			for (int i=0; i<str.length(); i++) {
				String tok = str.substring(i, i+1);  
				//System.out.println(tok); 
				// if (i==0) { c = cNum(tok); }
				if (i==1) { v = vNum(tok); }
				else if (i==2) { c = c2Num(tok); }
				//System.out.println(tok); 
			}
			//System.out.println(v + " " + c); 
			int count = ar[v][c]; 
			count++; 
			ar[v][c] = count; 
			c = -1; v = -1; 

			str = d.readLine();
		}

		/* counting VC */ 
		Hashtable<String, Integer> hashVC = new Hashtable<String, Integer>();

		sum = 0;
		sumsum = 0; 

		for (int i=0; i<nucleus; i++) { // # of V;
			sum = 0; 
			for (int j=0; j<coda; j++) { // # of C; 
				sum +=  ar[i][j]; 
				sumsum += ar[i][j]; 
			}
			String iString = vString(i); 
			hashVC.put(iString, sum); 
			//System.out.println(iString + ":" + sum); 
		}
		sum = 0; 
		for (int j=0; j<coda; j++) { // # of C; 
			sum = 0; 
			for (int i=0; i<nucleus; i++) { // # of V;
				sum +=  ar[i][j]; 
			}
			String jString = c2String(j); 
			hashVC.put(jString, sum); 
			//System.out.println(jString + ":" + sum); 

		}

		//System.out.println(sumsum); 


		// printing  v c numbers; 
		// 0 0 = # of ㅏ ㄱ 
		double rphi6VCsum = 0;   
		double ficherVCsum = 0; 
		int rphi6VCnum = 0;

		double[] ficherArrayVC = new double[nucleus*coda];
		double[] numArrayVC = new double[nucleus*coda]; 

		for (int i=0; i<nucleus; i++) { // # of V;
			for (int j=0; j<coda; j++) { // # of C; 
				String iString = vString(i);
				String jString = c2String(j);

				int anum = ar[i][j]; 
				int bnum = hashVC.get(iString) - anum; 
				int cnum = hashVC.get(jString) - anum; 
				int dnum = sumsum - (anum+bnum+cnum);


				double tpp = (double) anum / (anum + cnum);
				double deltapp = ((double) anum / (anum + cnum)) - ( (double) bnum / (bnum + dnum));
				double deltap = ((double) anum / (anum + bnum)) - ( (double) cnum / (cnum + dnum));

				double rphi5 = (double) ((double) (anum*dnum) - (bnum*cnum)) / Math.sqrt((anum + bnum) * (cnum + dnum) * (anum + cnum) * (bnum + dnum));
				double rphi6 = (double) Math.sqrt( deltap * deltapp);


                                FisherExact fe = new FisherExact(anum+bnum+cnum+dnum);
                                double ficherVCdouble = fe.getCumlativeP(anum,bnum,cnum,dnum);


				if (rphi6 < 0) rphi6 = rphi6 * -1;
				System.out.println( iString + "" + jString + "\t"  + rphi6 + "\t" + ficherVCdouble); 

				if (!Double.isNaN(rphi6)) numArrayVC[rphi6VCnum] =  rphi6;
				else numArrayVC[rphi6VCnum] =  0;

				if (!Double.isNaN(ficherVCdouble)) ficherArrayVC[rphi6VCnum] =  ficherVCdouble;
				else ficherArrayVC[rphi6VCnum] =  0;


				rphi6VCnum++; 
				if (!Double.isNaN(rphi6)) rphi6VCsum += rphi6; 	
				if (!Double.isNaN(ficherVCdouble)) ficherVCsum += ficherVCdouble;


				//System.out.println( iString + "" + jString + "\t:" + "\t" + anum + "\t" + bnum + "\t" + cnum + "\t" + dnum ); 
			}
		}

		System.out.println(); 
                System.out.println();


		d.close();
		double cvAverage = (double) rphi6CVsum / rphi6CVnum; 
		double vcAverage = (double) rphi6VCsum / rphi6VCnum; 

		System.out.println(fileName +",av\t" + "cv :" + cvAverage + "\t" + "vc :" + vcAverage); 
		double cvSD = calculateSD(numArrayCV); 
		double vcSD = calculateSD(numArrayVC); 
		System.out.println(fileName +",sd\t" + "cv :" + cvSD + "\t" + "vc :" + vcSD);


		cvAverage = (double) ficherCVsum / rphi6CVnum;
		vcAverage = (double) ficherVCsum / rphi6VCnum;
		System.out.println("ficher,av\t" + "cv :" + cvAverage + "\t" + "vc :" + vcAverage);

		cvSD = calculateSD(ficherArrayCV); 
		vcSD = calculateSD(ficherArrayVC);
		System.out.println("fichier,sd\t" + "cv :" + cvSD + "\t" + "vc :" + vcSD);

	}


	public static int cNum (String c) {
		int cnum = -1; 
		// onset = 18
		if (c.equals("ㄱ")) cnum = 0; 
		else if (c.equals("ㄴ")) cnum = 1;
		else if (c.equals("ㄷ")) cnum = 2;
		else if (c.equals("ㄹ")) cnum = 3;
		else if (c.equals("ㅁ")) cnum = 4;
		else if (c.equals("ㅂ")) cnum = 5;
		else if (c.equals("ㅅ")) cnum = 6;
		else if (c.equals("ㅈ")) cnum = 7;
		else if (c.equals("ㅊ")) cnum = 8;
		else if (c.equals("ㅋ")) cnum = 9;
		else if (c.equals("ㅌ")) cnum = 10;
		else if (c.equals("ㅍ")) cnum = 11;
		else if (c.equals("ㅎ")) cnum = 12;
		else if (c.equals("ㄲ")) cnum = 13;
		else if (c.equals("ㄸ")) cnum = 14;
		else if (c.equals("ㅃ")) cnum = 15;
		else if (c.equals("ㅆ")) cnum = 16;
		else if (c.equals("ㅉ")) cnum = 17;

		return cnum; 
	}

	public static int c2Num (String c) {
		int cnum = -1; 
		//coda = 7
		if (c.equals("ㄱ")) cnum = 0; 
		else if (c.equals("ㄴ")) cnum = 1;
		else if (c.equals("ㄷ")) cnum = 2;
		else if (c.equals("ㄹ")) cnum = 3;
		else if (c.equals("ㅁ")) cnum = 4;
		else if (c.equals("ㅂ")) cnum = 5;
		else if (c.equals("ㅅ")) cnum = 2; // ㄷ
		else if (c.equals("ㅇ")) cnum = 6;
		else if (c.equals("ㅈ")) cnum = 2; // ㄷ
		else if (c.equals("ㅊ")) cnum = 2; // ㄷ
		else if (c.equals("ㅋ")) cnum = 0; // ㄱ
		else if (c.equals("ㅌ")) cnum = 2; // ㄷ
		else if (c.equals("ㅍ")) cnum = 5; // ㅂ
		else if (c.equals("ㅎ")) cnum = 2; // ㄷ 
		else if (c.equals("ㄲ")) cnum = 0; // ㄱ
		else if (c.equals("ㄸ")) cnum = 2; // ㄷ
		else if (c.equals("ㅃ")) cnum = 5; // ㅂ
		else if (c.equals("ㅆ")) cnum = 2; // ㅅ
		else if (c.equals("ㅉ")) cnum = 2; // ㄷ
		//1) ㄳ : ㄱ
		else if (c.equals("ㄳ")) cnum = 0; // ㄱ
		//2) ㄵ : ㄴ
		else if (c.equals("ㄵ")) cnum = 1; // ㄴ
		//3) ㄽ : ㄹ
		else if (c.equals("ㄽ")) cnum = 3; // ㄹ
		//4) ㄾ: ㄹ
		else if (c.equals("ㄾ")) cnum = 3; // ㄹ
		//5) ㄺ: ㄱ ; 8) ㄺ: ㄱ
		else if (c.equals("ㄺ")) cnum = 0; // ㄱ
		//6) ㄼ: ㄹ  
		else if (c.equals("ㄼ")) cnum = 3; // ㄹ
		//7) ㅄ: ㅂ
		else if (c.equals("ㅄ")) cnum = 5; // ㅂ
		//9) ㄻ: ㅁ
		else if (c.equals("ㄻ")) cnum = 4; // ㅁ
		//10) ㄿ: ㅂ
		else if (c.equals("ㄿ")) cnum = 5; // ㅂ
		//11) ㄶ: ㄴ
		else if (c.equals("ㄶ")) cnum = 1; // ㄴ
		//12) ㅀ: ㄹ		
		else if (c.equals("ㅀ")) cnum = 3; // ㄹ

		else if (c.equals("ㄲ")) cnum = 0; // ㄱ
		else if (c.equals("ㄸ")) cnum = 2; // ㄷ
		else if (c.equals("ㅃ")) cnum = 5; // ㅂ
		else if (c.equals("ㅆ")) cnum = 2; // ㄷ
		else if (c.equals("ㅉ")) cnum = 2; // ㄷ



		return cnum; 
	}

	public static String c2String (int c) {
		String cnum = "O";  
		if (c==0) cnum = "ㄱ";
		else if (c==1) cnum = "ㄴ";
		else if (c==2) cnum = "ㄷ";
		else if (c==3) cnum = "ㄹ";
		else if (c==4) cnum = "ㅁ";
		else if (c==5) cnum = "ㅂ";
		else if (c==6) cnum = "ㅇ";

		return cnum;
	}


	public static String cString (int c) {
		String cnum = "O";  
		if (c==0) cnum = "ㄱ";
		else if (c==1) cnum = "ㄴ";
		else if (c==2) cnum = "ㄷ";
		else if (c==3) cnum = "ㄹ";
		else if (c==4) cnum = "ㅁ";
		else if (c==5) cnum = "ㅂ";
		else if (c==6) cnum = "ㅅ";
		else if (c==7) cnum = "ㅈ";
		else if (c==8) cnum = "ㅊ";
		else if (c==9) cnum = "ㅋ";
		else if (c==10) cnum = "ㅌ";
		else if (c==11) cnum = "ㅍ";
		else if (c==12) cnum = "ㅎ";
		else if (c==13) cnum = "ㄲ";
		else if (c==14) cnum = "ㄸ";
		else if (c==15) cnum = "ㅃ";
		else if (c==16) cnum = "ㅆ";
		else if (c==17) cnum = "ㅉ";

		return cnum;
	}

	public static String vString (int v) {
		String vnum = "O"; 
		if (v==0) vnum = "ㅏ"; 
		else if (v==1) vnum = "ㅓ";
		else if (v==2) vnum = "ㅗ"; 
		else if (v==3) vnum = "ㅜ"; 
		else if (v==4) vnum = "ㅡ"; 
		else if (v==5) vnum = "ㅣ"; 
		else if (v==6) vnum = "ㅔ/ㅐ"; 
		else if (v==7) vnum = "ㅢ";
		else if (v==8) vnum = "ㅝ";
		else if (v==9) vnum = "ㅘ";
		else if (v==10) vnum = "ㅚ/ㅞ/ㅙ";
		else if (v==11) vnum = "ㅟ";
		else if (v==12) vnum = "ㅖ/ㅒ";
		else if (v==13) vnum = "ㅑ";
		else if (v==14) vnum = "ㅛ";
		else if (v==15) vnum = "ㅠ";
		else if (v==16) vnum = "ㅕ";

		//else if (v==17) vnum = "ㅐ";
		//else if (v==18) vnum = "ㅞ";
		//else if (v==19) vnum = "ㅙ";
		//else if (v==20) vnum = "ㅒ";

		return vnum; 
	}

	public static int vNum (String v) {
		int vnum = -1; 
		if (v.equals("ㅏ")) vnum = 0; 
		else if (v.equals("ㅓ")) vnum = 1;
		else if (v.equals("ㅗ")) vnum = 2; 
		else if (v.equals("ㅜ")) vnum = 3; 
		else if (v.equals("ㅡ")) vnum = 4; 
		else if (v.equals("ㅣ")) vnum = 5; 
		else if (v.equals("ㅔ")) vnum = 6; 
		else if (v.equals("ㅐ")) vnum = 6; //ㅔ 
		else if (v.equals("ㅢ")) vnum = 7;
		else if (v.equals("ㅝ")) vnum = 8;
		else if (v.equals("ㅘ")) vnum = 9;
		else if (v.equals("ㅚ")) vnum = 10;
		else if (v.equals("ㅞ")) vnum = 10; //ㅚ
		else if (v.equals("ㅙ")) vnum = 10; //ㅚ
		else if (v.equals("ㅟ")) vnum = 11;
		else if (v.equals("ㅖ")) vnum = 12;
		else if (v.equals("ㅒ")) vnum = 12; // ㅖ
		else if (v.equals("ㅑ")) vnum = 13;
		else if (v.equals("ㅛ")) vnum = 14;
		else if (v.equals("ㅠ")) vnum = 15;
		else if (v.equals("ㅕ")) vnum = 16;

		return vnum; 
	}
	public static double calculateSD(double numArray[])
	{
		double sum = 0.0, standardDeviation = 0.0;
		int length = numArray.length;

		for(double num : numArray) {
			sum += num;
		}

		double mean = sum/length;

		for(double num: numArray) {
			standardDeviation += Math.pow(num - mean, 2);
		}

		return Math.sqrt(standardDeviation/length);
	}

	public static int factorial(int n){    
		if (n == 0)    
			return 1;    
		else    
			return(n * factorial(n-1));    
	}
	/*
	   BigDecimal fact = BigDecimal.valueOf(1);
	   for (int i = 1; i <= 8785856; i++)
	   fact = fact.multiply(BigDecimal.valueOf(i));
	   System.out.println(fact);
	 */
}
