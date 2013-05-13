package unitTests;

import org.junit.Test;

import qft.CoeffAdaptor;

import uk.ac.cam.cal56.maths.Complex;

public class CoeffAdaptorTest {

    @Test
    public void test() {

        // Tests the static function
        
        Complex[] passIn = new Complex[10];
        for(int i=0; i<10; i++) {
            passIn[i] = Complex.zero();
        }
        passIn[6] = Complex.one();
        
        Complex[] out = new Complex[21];
        
        out = CoeffAdaptor.setCoeffs(3, 21, passIn);
        
        assert(out[0].isZero() == true);
        assert(out[1].isZero() == true);
        assert(out[2].isZero() == true);
        assert(out[3].isZero() == true);
        assert(out[4].isZero() == true);
        assert(out[5].isZero() == true);
        assert(out[6].isZero() == true);
        assert(out[7].isZero() == true);
        assert(out[8].equals(new Complex(1.0, 0.0)));
        assert(out[9].isZero() == true);
        assert(out[10].isZero() == true);
        assert(out[11].isZero() == true);
        assert(out[12].isZero() == true);
        assert(out[13].isZero() == true);
        assert(out[14].isZero() == true);
        assert(out[15].isZero() == true);
        assert(out[16].isZero() == true);
        assert(out[17].isZero() == true);
        assert(out[18].isZero() == true);
        assert(out[19].isZero() == true);
        assert(out[20].isZero() == true);
        
    }

}
