import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

public class An_et_al_main {
	private static double alpha;
	private static double[][] fixedCharge;
	private static double[][] q;
	private static int nVar;
	private static double[][] flows;
	private static int P;
	private static double[][] distances;
	private static double Qmax;
	private static double rho = 1;

	/** Main
	 * 
	 * @param arg
	 * @throws FileNotFoundException
	 * @throws InterruptedException 
	 */
	public static void main(String[] arg) throws FileNotFoundException, InterruptedException{
		
		int[] N = {10,15};
		int[] hub = {3,5,7};
		double[] discount = {0.2};
		String[] failure = {
				"01-05",
				"05-10",
				"10-15",
				"15-20",
				"20-25",
				"01-10",
				"01-15",
				"01-20",
				"01-25"				
		};
		
		for(int n : N){
			for (int h : hub){
				for (double d : discount){
					for (String s:failure){
//						run(n, h, d, s);	
					}
				}
			}
		}
		
		run(15, 5, 0.2, "01-20");
	}

	/**
	 * 
	 * @param N
	 * @param hubs
	 * @param discount
	 * @param failure
	 * @throws InterruptedException
	 * @throws FileNotFoundException
	 */
	private static void run(int N, int hubs, double discount, String failure ) throws InterruptedException, FileNotFoundException{
		alpha = discount;
		fixedCharge = MyArray.read("Datasets/CAB/CAB" + N + "/fixedcharge.txt");
		q = MyArray.read("Datasets/CAB/Failures/" + failure + "/failures.txt");
		nVar = fixedCharge.length;
		flows = MyArray.read("Datasets/CAB/CAB" + N + "/Flows.txt");
		distances = MyArray.read("Datasets/CAB/CAB" + N + "/Distances.txt");
		P = hubs; // number of hubs to be located
		
		Qmax = getQmax();
		
		/** Filling in the flows matrix assymetrically */
//		for (int i = 0; i < nVar; i++) {
//			for (int j = 0; j < nVar; j++) {
//				flows[i][j] = tmpFlows[i][j] + tmpFlows[j][i];
//			}
//		}
		
		PrintWriter output = new PrintWriter(new File("D:/model_" + N + "_" + hubs + "_" + discount + "_" + failure + ".lp"));
				
		/**
		 * Objective Functions
		 */
		output.append("Minimize\n");
		
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int k=0;k<nVar;k++){
					for (int m=0;m<nVar;m++){
						if (k!=i && j!=m){
							double coef = Fikmj(i, k, m, j)*flows[i][j]*(1-q[k][0]-q(k, m));
							output.append(" + " + coef + " x" + i +"_" + k + "_" + m + "_" + j + " \n");
						}
					}
				}
			}
		}

		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				
				for (int m=0;m<nVar;m++){
					if (m!=j){
						double coef = Fikmj(i, i, m, j) * flows[i][j] * (1-q(i, m));
						output.append(" + " + coef + " x" + i +"_" + i + "_" + m + "_" + j + " \n");
					}
				}
				
				for (int k=0;k<nVar;k++){
					if (k!=i){
						double coef = Fikmj(i, k, j, j) * flows[i][j] * (1-q(j, k));
						output.append(" + " + coef + " x" + i +"_" + k + "_" + j + "_" + j + " \n");
					}
				}
				
				double coef = Fikmj(i, i, j, j) * flows[i][j];
				output.append(" + " + coef + " x" + i +"_" + i + "_" + j + "_" + j + " \n");
			}
		}
		
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int n=0;n<nVar;n++){
					output.append(" + omega"+i+"_"+j+"_"+ n + " - " + sigma(i,j) + " U"+i+"_"+j+"_"+n+" \n");
				}				
			}
		}
		
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int n=0;n<nVar;n++){
					output.append(" + theta"+i+"_"+j+"_"+ n + " - " + sigma(i,j) + " V" + i + "_" + j + "_" + n + " \n");
				}				
			}
		}

		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int n=0;n<nVar;n++){
					output.append(" + gamma"+i+"_"+j+"_"+ n + " - " + sigma(i,j)+" U"+i+"_"+j+"_"+n+" \n");
				}				
			}
		}
		
		/*for (int i = 0 ; i <nVar ; i++){
			output.append(" + " + fixedCharge[i][0] + " y" + i + " \n");
		}*/
	
		output.append("Subject to\n");
		
		/**
		 * Constraint 19
		 */
		for (int i=0; i<nVar; i++){ 
			for (int j=i+1; j<nVar;j++){
				
				for (int k=0;k<nVar;k++){
					output.append(" + x"+i+"_"+k+"_"+j+"_"+j);
				}
				output.append(" - y"+j+" = 0\n");
				
			}
		}
		
		/**
		 * Constraint 20
		 */
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar; j++){
				
				for (int m=0; m<nVar; m++){
					output.append(" + x" + i +"_"+ i +"_" + m + "_" + j);
				}
				output.append(" - y"+i+" = 0\n");				
			}
		}
		
		/**
		 * Constraint 21
		 */
		for (int i=0;i<nVar;i++){
			output.append(" + y"+i);
		}
		output.append(" = "+P+" \n");
		
		/**
		 * Constraint 22
		 */
		
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				
				for (int k=0;k<nVar;k++){
					for (int m=0;m<nVar;m++){
						output.append(" + x"+i+"_"+k+"_"+m+"_"+j);
					}
				}
				output.append(" = 1\n");
			}
		}
		
		/**
		 * Constraint 23
		 */
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int k=0;k<nVar;k++){
					
					output.append(" + U" + i + "_" + j + "_" + k);
					for (int m=0;m<nVar;m++){
						output.append(" + x" + i + "_" + k + "_" + m + "_" + j);
					}
					output.append(" - y" + k + " <= 0\n");
					
				}
			}
		}
		
		/**
		 * Constraint 24
		 */
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				
				for (int k=0; k<nVar; k++){
					output.append(" + U" + i + "_" + j + "_" + k);
				}
				
				for (int m=0;m<nVar;m++){
					output.append(" + x"+i+"_"+i+"_"+m+"_"+j);
					output.append(" + x"+i+"_"+j+"_"+m+"_"+j);
				}				
				output.append(" = 1\n");
				
			}
		}
		
		/**
		 * Constraint 25
		 */
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int m=0;m<nVar;m++){
					
					output.append(" + V"+i+"_"+j+"_"+m);
					for (int k=0;k<nVar;k++){
						output.append(" + x"+i+"_"+k+"_"+m+"_"+j);
					}
					output.append(" - y"+m+" <= 0\n");
				}
			}
		}
		
		/**
		 * Constraint 26
		 */
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				
				for (int m=0;m<nVar;m++){
					output.append(" + V"+i+"_"+j+"_"+m);
				}
				
				for (int k=0;k<nVar;k++){
					output.append(" + x"+i+"_"+k+"_"+j+"_"+j);
					output.append(" + x"+i+"_"+k+"_"+i+"_"+j);
				}
				
				output.append(" = 1\n");
				
			}
		}
		
		/**
		 * Constraint A1
		 */
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int n=0;n<nVar;n++){
					
					for (int k=0;k<nVar;k++){
						for (int m=0;m<nVar;m++){
							if (m!=k){
								double coef = rho * flows[i][j] * q[k][0]	* Fikmj(i, n, m, j);
								output.append(" + " + coef + " x"+i+"_"+k+"_"+m+"_"+j);
							}
						}
					}
					output.append(" - s"+i+"_"+j+"_"+n);
					output.append(" - omega"+i+"_"+j+"_"+n);
					output.append(" = - " + sigma(i,j) + " \n");
				}
			}
		}
		
		/**
		 * Constraint A2
		 */
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int n=0;n<nVar;n++){
					double coef = miu(i,j) + sigma(i,j);
					output.append(" + s"+i+"_"+j+"_"+n);
					output.append(" + "+ coef+ " U"+i+"_"+j+"_"+n);
					output.append(" <= "+coef+ " \n");
				}
			}
		}
		
		/**
		 * Constraint A3
		 */
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int n=0;n<nVar;n++){
					
					for (int k=0;k<nVar;k++){
						for (int m=0;m<nVar;m++){
							if (m!=k){
								double coef = rho * flows[i][j] * q[m][0]	* Fikmj(i, k, n, j);
								output.append(" + " + coef + " x"+i+"_"+k+"_"+m+"_"+j);
							}
						}
					}
					output.append(" - t"+i+"_"+j+"_"+n);
					output.append(" - theta"+i+"_"+j+"_"+n);
					output.append(" = - "+sigma(i,j)+ " \n");
				}
			}
		}
		
		/**
		 * Constraint A4
		 */
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int n=0;n<nVar;n++){
					double coef = miu(i,j) + sigma(i,j);
					output.append(" + t"+i+"_"+j+"_"+n);
					output.append(" + "+ coef+ " V"+i+"_"+j+"_"+n);
					output.append(" <= "+coef+ " \n");
				}
			}
		}
		
		/**
		 * Constraint A5
		 */
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int n=0;n<nVar;n++){
				
				for (int k=0;k<nVar;k++){
					double coef = rho * flows[i][j] * q[k][0] * Fikmj(i, n, n, j);
					output.append(" + " + coef + " x" + i + "_" + k + "_" + k + "_" + j);
				}
				output.append(" - r"+ i +"_"+ j +"_"+ n );
				output.append(" - gamma"+i+"_"+j+"_"+ n );
				output.append(" = - "+sigma(i,j)+ " \n");
				}
			}
		}
			
		/**
		 * Constraint A6
		 */
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int n=0;n<nVar;n++){
					double coef = miu(i,j) + sigma(i,j);
					output.append(" + r"+i+"_"+j+"_"+n);
					output.append(" + "+ coef+ " U"+i+"_"+j+"_"+n);
					output.append(" <= "+coef+ " \n");
				}
			}
		}
		
		/**
		 * Binaries
		 */
		output.append("Binaries\n");
		
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int k=0;k<nVar;k++){
					for (int m=0;m<nVar;m++){
						output.append("x" + i + "_" + k + "_" + m + "_" + j + " \n");
					}
				}
			}
		}
		
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int n=0 ; n<nVar ;n++){
					output.append("U" + i + "_" + j + "_" + n + " \n");
					output.append("V" + i + "_" + j + "_" + n + " \n");
				}
			}
		}
		
		for (int i=0;i<nVar;i++){
			 output.append("y" + i + " \n");
		}
		
		
		output.close();
		System.out.println("Model built.");
		
		ProcessBuilder pb = new ProcessBuilder("cmd", "/c", "start"
				 , "gurobi_cl"
				 , "ResultFile=D:/Results_" + N + "_" + hubs + "_" + discount + "_" + failure + ".sol"
				 , "D:/model_" + N + "_" + hubs + "_" + discount + "_" + failure + ".lp");
		
		try {
			 pb.start().waitFor();
			 } catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	/** qkm */
	private static double q(int k, int m) {
		if (m==k)
			return 0;
		else 
			return q[m][0];
	}
	
	/** Fikmj */
	private static double Fikmj(int i, int k, int m, int j) {
		double cost = distances[i][k] + (1 - alpha) * distances[k][m]
				+ distances[m][j];
		return cost;
	}
	
	/** Sigma */
	private static double sigma(int i, int j){
		return 0;
	}
	
	/** miu */
	private static double miu(int i, int j){
		double output = rho * flows[i][j] * Fmax(i,j) * Qmax;
		return output;
	}

	/** Fmax
	 * 
	 * @param i
	 * @param j
	 * @return
	 */
	private static double Fmax(int i, int j){
		double y=0;
		for (int k=0;k<nVar;k++){
			for (int m=0;m<nVar;m++){
				double cost = Fikmj(i, k, m, j);
				if (y<cost)
					y=cost;
			}
		}
		return y;
	}
	
	private static double getQmax (){
		double output = 0;
		for (int i = 0 ; i < q.length ; i++)
			output = ( q[i][0] > output ) ? q[i][0] : output;
		return output;
	}
}
