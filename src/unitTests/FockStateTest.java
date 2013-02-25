package unitTests;

import static org.junit.Assert.*;

import org.junit.Test;

import qft.FockState;

public class FockStateTest {
    
    @Test
    public void testSetUpAndGet() {
        
        //Test constructor, setting as vacuum, setters and getters
        
        FockState test1 = new FockState(10);
        double epsilon = 0.0000001;
        assertEquals(test1.getCoeff(7), 0, epsilon);
        assertEquals(test1.getCoeff(0), 0, epsilon);
        assertEquals(test1.getCoeff(10), -1, epsilon);
        test1.setCoeff(3, 5);
        test1.setCoeff(10, 3);
        assertEquals(test1.getCoeff(10), -1, epsilon);
        assertEquals(test1.getCoeff(3), 5, epsilon);
        test1.setAsVacuum();
        assertEquals(test1.getCoeff(3), 0, epsilon);
        assertEquals(test1.getSystemSize(), 10, epsilon);
        
        
    }

    @Test
    public void testCreatAndAnnih() {
        
        //Tests creatin and annihilation operators (and validity checks)
        
        FockState test2 = new FockState(13);
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
    public void testIncrement() {
        
        //Tests the state incrementation method in various cases
        
        FockState test3 = new FockState(5);
        
        test3.setCoeff(0, 5);
        test3.incrementState();
        assertEquals("simple check 1", test3.getCoeff(0), 4);
        assertEquals("simple check 2", test3.getCoeff(1), 1);
        test3.incrementState();
        assertEquals("simple check 3", test3.getCoeff(0), 4);
        assertEquals("simple check 4", test3.getCoeff(2), 1);
        
        test3.setAsVacuum();
        test3.setCoeff(2, 2);
        test3.setCoeff(4, 2);
        test3.incrementState();
        assertEquals("middle check 1", test3.getCoeff(2), 1);
        assertEquals("middle check 2", test3.getCoeff(3), 3);
        
        test3.setAsVacuum();
        test3.setCoeff(3, 2);
        test3.setCoeff(4, 1);
        test3.incrementState();
        assertEquals("end check 1", test3.getCoeff(3), 1);
        assertEquals("end check 2", test3.getCoeff(4), 2);
        
        test3.setAsVacuum();
        test3.setCoeff(4, 2);
        test3.incrementState();
        assertEquals("particle increase check 1", test3.getCoeff(0), 3);
        assertEquals("particle increase check 2", test3.getCoeff(4), 0);
        
    }

    @Test
    public void testIndexing() {
        
        //Tests the setting of a state according to index and calculating index according to state
        
        double epsilon = 0.000001;
        FockState test4 = new FockState(2);
        FockState test5 = new FockState(4);
        FockState test6 = new FockState(7);
        
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
        
        FockState test7 = new FockState(65);
        FockState test8 = new FockState(65);
        double epsilon = 0.000001;
        test7.setCoeff(53, 3);
        test8.makeSameAs(test7);
        assertEquals(test7.getCoeff(53), test8.getCoeff(53), epsilon);
        assertEquals(test7.getCoeff(3), test8.getCoeff(3), epsilon);
        
    }
    
}
