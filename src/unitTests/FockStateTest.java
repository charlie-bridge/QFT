package unitTests;

import static org.junit.Assert.*;

import org.junit.Test;

import qft.FockState;

public class FockStateTest {
    
    @Test
    public void testSetUpAndGet() {
        
        //Test constructor, setting as vacuum, setters and getters
        
        FockState test1 = new FockState(10, 0.1, 1.0);
        double epsilon = 0.0000001;
        assertEquals(test1.getCoeff(7), 0, epsilon);
        assertEquals(test1.getCoeff(0), 0, epsilon);
        assert(test1.getCoeff(10) == null);
        test1.setCoeff(3, 5);
        test1.setCoeff(10, 3);
        assert(test1.getCoeff(10) == null);
        assertEquals(test1.getCoeff(3), 5, epsilon);
        test1.setAsVacuum();
        assertEquals(test1.getCoeff(3), 0, epsilon);
        assertEquals(test1.getSystemSize(), 10, epsilon);
        
        
    }

    @Test
    public void testCreateAndAnnihil() {
        
        //Tests creation and annihilation operators (and validity checks)
        
        FockState test2 = new FockState(13, 0.1, 1.0);
        double epsilon = 0.000001;
        test2.setCoeff(4, 2);
        test2.setCoeff(3, 1);
        test2.applyCreation(11);
        assertEquals(test2.getCoeff(11), 1, epsilon);
        test2.applyCreation(11);
        assertEquals(test2.getCoeff(11), 2, epsilon);
        assert(test2.isValid() == true);
        test2.applyAnnihil(4);
        test2.applyAnnihil(3);
        assertEquals(test2.getCoeff(4), 1, epsilon);
        assertEquals(test2.getCoeff(3), 0, epsilon);
        assert(test2.isValid() == true);
        test2.applyAnnihil(4);
        assertEquals(test2.getCoeff(4), 0, epsilon);
        assert(test2.isValid() == true);
        test2.applyAnnihil(3);
        assertEquals(test2.getCoeff(3), -1, epsilon);
        assert(test2.isValid() == false);
        test2.applyCreation(3);
        assertEquals(test2.getCoeff(3), 0, epsilon);
        assert(test2.isValid() == false);
        
    }
    
    @Test
    public void testIndexing() {
        
        //Tests the setting of a state according to index and calculating index according to state
        
        double epsilon = 0.000001;
        FockState test4 = new FockState(2, 0.1, 1.0);
        FockState test5 = new FockState(4, 0.1, 1.0);
        FockState test6 = new FockState(7, 0.1, 1.0);
        
        test4.setCoeff(0, 3);
        test4.setCoeff(1, 2);
        assertEquals(test4.calcIndex(), 17, epsilon);
        
        test4.setAsIndex(11);
        assertEquals(test4.getCoeff(0), 3, epsilon);
        assertEquals(test4.getCoeff(1), 1, epsilon);
        
        test5.setCoeff(0, 1);
        test5.setCoeff(2, 1);
        assertEquals(test5.calcIndex(), 8, epsilon);
        
        test5.setAsVacuum();
        test5.setCoeff(1, 1);
        test5.setCoeff(2, 1);
        assertEquals(test5.calcIndex(), 14, epsilon);
        
        test5.setAsIndex(14);
        assertEquals(test5.getCoeff(0), 0, epsilon);
        assertEquals(test5.getCoeff(1), 1, epsilon);
        assertEquals(test5.getCoeff(2), 1, epsilon);
        assertEquals(test5.getCoeff(3), 0, epsilon);
        
        test6.setCoeff(0, 6);
        assertEquals(test6.calcIndex(), 21, epsilon);
        
        test6.setAsVacuum();
        test6.setAsIndex(20);
        assertEquals(test6.getCoeff(0), 0, epsilon);
        assertEquals(test6.getCoeff(1), 0, epsilon);
        assertEquals(test6.getCoeff(2), 0, epsilon);
        assertEquals(test6.getCoeff(3), 1, epsilon);
        assertEquals(test6.getCoeff(4), 0, epsilon);
        assertEquals(test6.getCoeff(5), 0, epsilon);
        assertEquals(test6.getCoeff(6), 0, epsilon);
       
    }

    @Test
    public void testMakeSame() {
        
        //Tests the state comparison method
        
        FockState test7 = new FockState(65, 0.1, 1.0);
        FockState test8 = new FockState(65, 0.1, 1.0);
        double epsilon = 0.000001;
        test7.setCoeff(53, 3);
        test8.makeSameAs(test7);
        assertEquals(test7.getCoeff(53), test8.getCoeff(53), epsilon);
        assertEquals(test7.getCoeff(3), test8.getCoeff(3), epsilon);
        
    }
    
    @Test
    public void testEnergyCalc() {
    	
    	//Tests the energy calculation method
    	
    	FockState test9 = new FockState(6, 0.55, 2.2);
    	double epsilon = 0.0000001;
    	test9.setCoeff(2, 1);
    	test9.setCoeff(3, 2);
    	test9.setCoeff(4, 3);
    	assertEquals(test9.calcEnergy(), 22.0617339, epsilon);
    	
    }
    
    @Test
    public void testFreqCalc() {
    	
    	//Tests the static frequency calculation
    	//The frequency calculation for two consecutive indexes 
    	//(starting with an odd one) should be the same
    	
    	double test10 = FockState.calcFrequency(7, 204, 0.15, 0.7);
    	double epsilon = 0.000000001;
    	assertEquals(test10, 1.078764327, epsilon);
    	
    	assertEquals(FockState.calcFrequency(1, 6, 0.2, 1.0), FockState.calcFrequency(2, 6, 0.2, 1.0), epsilon);
    	
    }
   
    @Test
    public void testFourierGets() {

    	double epsilon = 0.00000001;
    	
    	Integer onePIndex = FockState.getIndex1PState(2, 4, 0.1, 2.0);
    	assertEquals(onePIndex, 20 ,epsilon);
    	
    	onePIndex = FockState.getIndex1PState(1, 4, 0.1, 2.0);
    	assertEquals(onePIndex, 5 ,epsilon);
    	
    	onePIndex = FockState.getIndex1PState(3, 4, 0.1, 2.0);
    	assertEquals(onePIndex, 2 ,epsilon);
    	
    	Integer twoPIndex = FockState.getIndex2PState(1, 1, 4, 0.1, 1.5);
    	assertEquals(twoPIndex, 54 ,epsilon);
    	
    	twoPIndex = FockState.getIndex2PState(0, 3, 4, 0.1, 1.5);
    	assertEquals(twoPIndex, 4 ,epsilon);
    	
    	twoPIndex = FockState.getIndex2PState(3, 0, 4, 0.1, 1.5);
    	assertEquals(twoPIndex, 4 ,epsilon);
    	
    }

}
