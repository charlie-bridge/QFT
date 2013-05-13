package unitTests;

import static org.junit.Assert.assertEquals;

import java.util.HashMap;
import java.util.Map;

import org.junit.Before;
import org.junit.Test;

import qft.GeneralState;
import uk.ac.cam.cal56.maths.Complex;
import uk.ac.cam.cal56.qft.Interaction;
import uk.ac.cam.cal56.qft.impl.MomentumWavePacket;

public class GeneralStateTest {

	private GeneralState _genState;
	Map<Interaction, Double> lambdas = new HashMap<Interaction, Double>();

	@Before
	public void setUp() throws Exception {
		
		lambdas.put(Interaction.PHI_CUBED, 0.5);
		_genState = new GeneralState(6, 3, 1.1, 0.1, 0.1, lambdas, new MomentumWavePacket(6));
		
	}

	@SuppressWarnings("unused")
	@Test
	public void test() {
	    
	    // Will it run test
		
		double epsilon = 0.0000001;
		_genState.step();
		_genState.step();
		_genState.step();
		_genState.step();
		_genState.step();
		_genState.step();
		assertEquals(_genState.getTime(),  0.6, epsilon);
		Complex zerop;
		zerop = _genState.getVacuum(); 
		Complex[] onep;
		onep = _genState.get1PMom();
		Complex[][] twop;
		twop = _genState.get2PMom();		
		
	}

}
