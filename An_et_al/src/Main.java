import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

public class Main {
	private static double alpha = 0.2;
	private static double[][] tmpFlows = MyArray.read("w.txt");
	private static int nVar = tmpFlows.length;
	private static double[][] flows = new double[nVar][nVar];
	// private static double[][] fixedCharge = MyArray.read("fixedcharge.txt");
	private static double[][] coordinates = MyArray.read("coordinates.txt");
	private static double[][] distances = Distance.get(coordinates);
	private static int P = 3; // number of hubs to be located
	private static double q = 0.05; // probability of a node being disrupted
//	private static int M = nVar * 4; // the big M
	private static double Qmax;
	private static double rho = 1;

	/** qkm */
	private static double q(int k, int m) {
		if (m==k)
			return 0;
		else 
			return q;
	}
	
	/** Fikmj */
	private static double Fijkm(int i, int k, int m, int j) {
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
				double cost = Fijkm(i, k, m, j);
				if (y<cost)
					y=cost;
			}
		}
		return y;
	}
	
	/** Main
	 * 
	 * @param arg
	 * @throws FileNotFoundException
	 */
	public static void main(String[] arg) throws FileNotFoundException{
		
		Qmax = q;
		
		/** Filling in the flows matrix assymetrically */
		for (int i = 0; i < nVar; i++) {
			for (int j = 0; j < nVar; j++) {
				flows[i][j] = tmpFlows[i][j] + tmpFlows[j][i];
			}
		}
		
		PrintWriter output = new PrintWriter(new File("C:/Users/Pouya/git/An_et_al/An_et_al/ModelAndResults/model.lp"));
				
		/**
		 * Objective Functions
		 */
		output.append("Minimize\n");
		
		for (int i=0; i<nVar; i++){
			for (int j=i+1; j<nVar;j++){
				for (int k=0;k<nVar;k++){
					for (int m=0;m<nVar;m++){
						if (k!=i && j!=m){
							double coef = Fijkm(i, k, m, j)*flows[i][j]*(1-q-q(k, m));
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
						double coef = Fijkm(i, i, m, j) * flows[i][j] * (1-q(i, m));
						output.append(" + " + coef + " x" + i +"_" + i + "_" + m + "_" + j + " \n");
					}
				}
				
				for (int k=0;k<nVar;k++){
					if (k!=i){
						double coef = Fijkm(i, k, j, j) * flows[i][j] * (1-q(j, k));
						output.append(" + " + coef + " x" + i +"_" + k + "_" + j + "_" + j + " \n");
					}
				}
				
				double coef = Fijkm(i, i, j, j) * flows[i][j];
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
					output.append(" + 2 gamma"+i+"_"+j+"_"+ n + " - " + sigma(i,j)+" U"+i+"_"+j+"_"+n+" \n");
				}				
			}
		}
	
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
								double coef = rho * flows[i][j] * q	* Fijkm(i, n, m, j);
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
								double coef = rho * flows[i][j] * q	* Fijkm(i, k, n, j);
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
					double coef = rho * flows[i][j] * q * Fijkm(i, n, n, j);
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
		System.out.println("done!");
		
		// Solve the model
		/*ProcessBuilder pb = new ProcessBuilder("cmd", "/c", "start",
				 "gurobi_cl",
				 "ResultFile = C:/Users/Pouya/git/An_et_al/An_et_al/ModelAndResults/Results222.sol"
				 ,"C:/Users/Pouya/git/An_et_al/An_et_al/ModelAndResults/model222.lp");
		
		ProcessBuilder pb = new ProcessBuilder("cmd", "/c", "start"
				 , "C:/gurobi603/win64/bin/gurobi_cl"
				 , "ResultFile=D:/Results.sol"
				 , "D:/model.lp");
		
		try {
			 pb.start();
			 } catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		
	}
}
