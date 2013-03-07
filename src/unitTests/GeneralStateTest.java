package unitTests;

import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

import qft.GeneralState;
import uk.ac.cam.cal56.maths.Complex;

public class GeneralStateTest {

	private GeneralState _genState;
	
	@Before
	public void setUp() throws Exception {
		
		_genState = new GeneralState(6, 3, 1.1, 0.1, 0.1, 1);
		
	}

	@SuppressWarnings("unused")
	@Test
	public void test() {
		
		double epsilon = 0.0000001;
		_genState.step();
		_genState.step();
		_genState.step();
		_genState.step();
		_genState.step();
		_genState.step();
		assertEquals(_genState.getTime(),  0.6, epsilon);
		Complex zerop;
		zerop = _genState.get0P(); 
		Complex[] onep;
		onep = _genState.get1PMom();
		Complex[][] twop;
		twop = _genState.get2PMom();
		
		
		
	}

}
