package unitTests;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import qft.FockState;
import qft.InteractionMatrix;

public class InteractionMatrixTest {
    
	@Test
	public void testIncrMomLabel() {
		
		//Tests the momenta label incrementation method
		
		InteractionMatrix test1 = new InteractionMatrix(10, 0.1, 1.2, 30, 3);
		InteractionMatrix test2 = new InteractionMatrix(50, 0.1, 1.2, 30, 5);
		int[] momenta1 = new int[2];
		int[] momenta2 = new int[4];
		double epsilon = 0.00000001;
		
		momenta1[0] = 2;
		momenta1[1] = 7;
		test1.incrementMomentaLabel(momenta1);
		assertEquals(momenta1[0], 3, epsilon);
		assertEquals(momenta1[1], 7, epsilon);
		momenta1[0] = 9;
		test1.incrementMomentaLabel(momenta1);
		assertEquals(momenta1[0], 0, epsilon);
		assertEquals(momenta1[1], 8, epsilon);
		test1.incrementMomentaLabel(momenta1);
		assertEquals(momenta1[0], 1, epsilon);
		assertEquals(momenta1[1], 8, epsilon);
		
		momenta2[0] = 49;
		momenta2[1] = 49;
		momenta2[2] = 49;
		momenta2[3] = 0;
		test2.incrementMomentaLabel(momenta2);
		assertEquals(momenta2[0], 0, epsilon);
		assertEquals(momenta2[1], 0, epsilon);
		assertEquals(momenta2[2], 0, epsilon);
		assertEquals(momenta2[3], 1, epsilon);
		
	}
    
	@Test
	public void testApplyOps() {
		
		//Tests the operator application method
		
		//TODO Maybe try some more, this is very important to a couple more test cases would be good
		
		double epsilon = 0.0000001;
		InteractionMatrix test3 = new InteractionMatrix(10, 0.1, 1.2, 50, 4);
		InteractionMatrix test4 = new InteractionMatrix(15, 0.2, 1.5, 40, 3);
		int[] momenta3 = new int[3];
		int[] momenta4 = new int[2];
		int opType3 =  15;
		int opType4 = 3;
		FockState testState3 = new FockState(10, 0.1, 1.2);
		FockState testState4 = new FockState(15, 0.2, 1.5);
		testState3.setCoeff(0, 2);
		testState3.setCoeff(1, 1);
		testState3.setCoeff(3, 3);
		testState3.setCoeff(6, 4);
		testState4.setCoeff(3, 1);
		testState4.setCoeff(4, 2);
		momenta3[0] = 3;
		momenta3[1] = 0;
		momenta3[2] = 1;
		momenta4[0] = 0;
		momenta4[1] = 1;
		
		
		test3.applyOperators(testState3, momenta3, opType3);
		test4.applyOperators(testState4, momenta4, opType4);
		
		assertEquals(testState3.getCoeff(0), 1, epsilon);
		assertEquals(testState3.getCoeff(1), 0, epsilon);
		assertEquals(testState3.getCoeff(2), 0, epsilon);
		assertEquals(testState3.getCoeff(3), 2, epsilon);
		assertEquals(testState3.getCoeff(4), 0, epsilon);
		assertEquals(testState3.getCoeff(5), 0, epsilon);
		assertEquals(testState3.getCoeff(6), 3, epsilon);
		assertEquals(testState3.getCoeff(7), 0, epsilon);
		assert(testState3.isValid() == true);
		assertEquals(test3.getFactors(), 4.89897948, epsilon);
	
		assertEquals(testState4.getCoeff(0), -1, epsilon);
		assertEquals(testState4.getCoeff(1), 0, epsilon);
		assertEquals(testState4.getCoeff(3), 1, epsilon);
		assertEquals(testState4.getCoeff(4), 2, epsilon);
		assert(testState4.isValid() == false);
		assertEquals(test4.getFactors(), 0.0, epsilon);
		
	}
	
	@Test
	public void testCalc() {
		
		//Tests the overall calculation method
		
		InteractionMatrix test6 = new InteractionMatrix(2, 0.11, 1.3, 3, 3);
		double epsilon = 0.0000000001;
		test6.calcMatrix();
		
		assertEquals(test6.getRow(0).get(1), 0.1362031092, epsilon);
		assertEquals(test6.getRow(0).get(6), 0.1038061339, epsilon);
		assertEquals(test6.getRow(1).get(0), 0.1362031092, epsilon);
		assertEquals(test6.getRow(1).get(3), 0.3724177824, epsilon);
	    assertEquals(test6.getRow(1).get(5), 0.01282278625, epsilon);
		assertEquals(test6.getRow(2).get(4), 0.1543372675, epsilon);
		assertEquals(test6.getRow(3).get(1), 0.3724177824, epsilon);
		assertEquals(test6.getRow(3).get(6), 0.6763228326, epsilon);
		assertEquals(test6.getRow(4).get(2), 0.1543372675, epsilon);
		assertEquals(test6.getRow(5).get(1), 0.01282278625, epsilon);
		assertEquals(test6.getRow(6).get(0), 0.1038061339, epsilon);
		assertEquals(test6.getRow(6).get(3), 0.6763228326, epsilon);
		
	}
	
}
