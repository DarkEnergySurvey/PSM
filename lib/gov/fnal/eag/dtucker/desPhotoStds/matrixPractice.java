/**
 * matrixPractice.java   The DES Collaboration  2006
 */
package gov.fnal.eag.dtucker.desPhotoStds;





import cern.colt.matrix.*;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.LUDecomposition;
import cern.colt.matrix.linalg.QRDecomposition;
import cern.colt.matrix.linalg.SingularValueDecomposition;


/**
 * @author dtucker
 *
 */
public class matrixPractice {

    public static void main (String[] args) {

    		DoubleMatrix2D matrix;
		matrix = new DenseDoubleMatrix2D(3,4);
		System.out.println(matrix); 
		int rows = matrix.rows();
		int columns = matrix.columns();
		System.out.println(rows + "\t" + columns);
		int row = 2, column = 0;
		matrix.set(row,column, 7);
		System.out.println(matrix.get(row,column));
		double sum = 0;
		for (row=rows; --row >= 0; ) {
		   for (column=columns; --column >= 0; ) {
		      sum += matrix.get(row,column); // bounds check
		      //sum += matrix.getQuick(row,column); // no bounds check
		   }
		}
		System.out.println(sum); 
		matrix.viewPart(1,0,2,2).assign(1);
		System.out.println(matrix); 
		System.out.println(matrix.viewDice());
		System.out.println(matrix.viewColumnFlip());
		System.out.println(matrix.viewRow(1));
		int[] rowIndexes = {0,2};
		int[] columnIndexes = {2,3,1,1};
		System.out.println(matrix.viewSelection(rowIndexes,columnIndexes));
		System.out.println(matrix.viewSorted(1).viewRowFlip());
		matrix.viewPart(0,1,2,2).viewRowFlip().viewColumn(0).assign(2);
		System.out.println(matrix); 
		
    		DoubleMatrix2D A = new DenseDoubleMatrix2D(3,3);
    		DoubleMatrix2D b = new DenseDoubleMatrix2D(3,1);
    		DoubleMatrix2D x = new DenseDoubleMatrix2D(3,1);
    		DoubleMatrix2D x1 = new DenseDoubleMatrix2D(3,1);
    		DoubleMatrix2D x2 = new DenseDoubleMatrix2D(3,1);
    		DoubleMatrix2D x3 = new DenseDoubleMatrix2D(3,1);
    		DoubleMatrix2D x4 = new DenseDoubleMatrix2D(3,1);
    		DoubleMatrix2D Ainv;
    		
    		double[][] array = {{1.,2.,3},{4.,5.,6.},{7.,8.,10.}};
    		for (int i=0; i<3; i++) {
    			for (int j=0; j<3; j++) {
    				A.set(i,j,array[i][j]);
    			}
    		}
    		System.out.println("A is \n" + A);
    		//for (int i=0; i<3; i++) {
    		//	b.set(i,0,Math.random());
    		//}
    		b.set(0,0,5.);
    		b.set(1,0,2.);
    		b.set(2,0,1.);
    		
    		System.out.println("b is \n" + b);
    		Algebra alg = new Algebra();
    		x = alg.solve(A,b);
    		System.out.println("x is \n" + x);
    	    System.out.println("Computing LU...");
    		LUDecomposition luq = new LUDecomposition(A);
    	    x1 = luq.solve(b);
    	    System.out.println("done!\nLU Decomposition, x:\n"+ x1);
  		QRDecomposition qr = new QRDecomposition(A);
    	    x2 = qr.solve(b);
    	    System.out.println("done!\nQR Decomposition, x:\n"+ x2);
    	    Ainv = alg.inverse(A);
    	    System.out.println("done!\nInverse of A is:\n" + Ainv);
    	    x3 = alg.mult(Ainv,b);
    	    System.out.println("done!\nAinv*b, x:\n"+ x3);
    	    System.out.println("Computing SVD...");
    	    
    	    DoubleMatrix2D U = new DenseDoubleMatrix2D(3,3);
    	    DoubleMatrix2D V = new DenseDoubleMatrix2D(3,3);
    	    DoubleMatrix2D S = new DenseDoubleMatrix2D(3,3);
    	    SingularValueDecomposition svd = new SingularValueDecomposition(A);
    	    System.out.println("done!\nRank of A:"+svd.rank());
    	    U = svd.getU();
    	    V = svd.getV();
    	    S = svd.getS();
    	    System.out.println("SVD U:\n" + U + "SVD V:\n" + V + "SVD S:\n" + S);
    	    System.out.println("SVD " + svd.toString());
  	
    	    
    	    
        DoubleMatrix2D AA = new DenseDoubleMatrix2D(2000,2000);
    		DoubleMatrix2D AAinv;
    		System.out.println("1");
    		for (int i=0; i<2000; i++) {
    			for (int j=0; j<2000; j++) {
    				AA.set(i,j,Math.random());
    			}
    		}
    		System.out.println("start inversion");
    	    AAinv = alg.inverse(AA);
    		System.out.println("finished");

    	    //System.out.println("\nAA:\n" + AA + "\nAAinv:\n" + AAinv);

    	    double[][] matrixArray = new double[10][10];
    	    for (int i=0; i<10; i++) {
    	    		for (int j=0; j<10; j++) {
    	    			matrixArray[i][j] = 0;
    	    		}
    	    }
    	    matrixArray[3][3] = 3;
    	    matrixArray[3][3] = matrixArray[3][3]+1;
    	    matrixArray[3][3] = matrixArray[3][3]+0.45;
    	    
    	    for (int i=0; i<10; i++) {
	    		for (int j=0; j<10; j++) {
	    			System.out.println(i + "\t" + j + "\t" + matrixArray[i][j]);
	    		}
    	    }
    }
    
}









