package unitTests;

import static org.junit.Assert.*;

import org.junit.Test;

import utilities.Complex;

public class ComplexTest {
    
    // COMPLEX CLASS SUPERCEDED BY CODEBASE VERSION, TESTS HERE ARE REDUNDANT

	@Test
	public void testSetUpAndGet() {
		
		Complex tester1 = new Complex(1.0, 2.0);
		double epsilon = 0.0000001;
		assertEquals(tester1.realPart(), 1.0, epsilon);
		assertEquals(tester1.complexPart(), 2.0, epsilon);
		
	}

	@Test
	public void testReset() {
		
		Complex tester2 = new Complex(3.0, 4.0);
		tester2.setReal(5.0);
		tester2.setComp(6.0);
		double epsilon = 0.0000001;
		assertEquals(tester2.realPart(), 5.0, epsilon);
		assertEquals(tester2.complexPart(), 6.0, epsilon);
		
	}
	
	@Test
	public void testAddAndSubtract() {
		
		Complex tester3 = new Complex(3.5, 7.2);
		Complex tester4 = new Complex(2.9, 8.6);
		double epsilon = 0.0000001;
		tester3.subtract(tester4);
		assertEquals(tester3.realPart(), 0.6, epsilon);
		assertEquals(tester3.complexPart(), -1.4, epsilon);
		tester3.add(tester4);
		assertEquals(tester3.realPart(), 3.5, epsilon);
		assertEquals(tester3.complexPart(), 7.2, epsilon);
		
	}

	@Test
	public void testModulus() {
		
		Complex tester5 = new Complex(3.0, 4.0);
		assert(tester5.modulus() == 5.0);
		
	}
	
}
