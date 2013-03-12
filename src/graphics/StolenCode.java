package graphics;

import java.awt.EventQueue;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import qft.GeneralState;
import uk.ac.cam.cal56.graphics.QFTSandbox;
import uk.ac.cam.cal56.qft.interactingtheory.Interaction;
import uk.ac.cam.cal56.qft.interactingtheory.impl.MomentumWavePacket;

@SuppressWarnings("serial")
public class StolenCode extends QFTSandbox {

	// Stolen QFTSandbox
	
	protected static final String FRAME_TITLE = "BLEEERURUHJELAK";
	
	private static final int _N_MAX = 256;
	
	
	@Override
	protected void setupSliders() {
	    
	    _NSlider = new JSlider(N_MIN, _N_MAX, N_DEFAULT);
	    
        // add calculate sliders
        setupGeneralSlider(_NSlider, _NValue, N_MIN, _N_MAX, 2, int.class, "Number of lattice points");
        setupGeneralSlider(_PmaxSlider, _PmaxValue, PMAX_MIN, PMAX_MAX, 3, int.class, "Number of particles considered");
        setupGeneralSlider(_dxSlider, _dxValue, encode(DX_MIN), encode(DX_MAX), 4, double.class,
                           "Lattice point separation");
        setupGeneralSlider(_mSlider, _mValue, encode(M_MIN), encode(M_MAX), 5, double.class, "Particle mass");
        // add real time sliders
        setupGeneralSlider(_dtSlider, _dtValue, encode(DT_MIN), encode(DT_MAX), 7, double.class, "Time step");
        setupGeneralSlider(_stepsSlider, _stepsValue, STEPS_MIN, STEPS_MAX, 8, int.class, "Steps calculated per frame");
        setupGeneralSlider(_lambdaSquaredSlider, _lambdaSquaredValue, encode(LAMBDA_MIN), encode(LAMBDA_MAX), 9,
                           double.class, "2-vertex interaction strength");
        setupGeneralSlider(_lambdaCubedSlider, _lambdaCubedValue, encode(LAMBDA_MIN), encode(LAMBDA_MAX), 10,
                           double.class, "3-vertex interaction strength");

        // then add real time update listeners (time step and interaction strength)
        _dtSlider.addChangeListener(new ChangeListener() { // update time step
            public void stateChanged(ChangeEvent e) {
                if (_state != null)
                    _state.setTimeStep(decode(_dtSlider.getValue()));
            }
        });
        _lambdaSquaredSlider.addChangeListener(new ChangeListener() { // update interaction strength
            public void stateChanged(ChangeEvent e) {
                if (_state != null)
                    _state.setInteractionStrength(Interaction.PHI_SQUARED, decode(_lambdaSquaredSlider.getValue()));
            }
        });
        _lambdaCubedSlider.addChangeListener(new ChangeListener() { // update interaction strength
            public void stateChanged(ChangeEvent e) {
                if (_state != null)
                    _state.setInteractionStrength(Interaction.PHI_CUBED, decode(_lambdaCubedSlider.getValue()));
            }
        });
    }
	
	@Override
	protected void setupQuantumState() {
		Map<Interaction, Double> lambdas = new HashMap<Interaction, Double>();
		
		//lambdas.put(Interaction.PHI_CUBED, 0.0);
		lambdas.put(Interaction.PHI_CUBED, decode(_lambdaCubedSlider.getValue()));

        int N = _NSlider.getValue();
        if (_wavepacket == null)
            _wavepacket = new MomentumWavePacket(N);
        _state = new GeneralState(N, _PmaxSlider.getValue(), decode(_mSlider.getValue()),
            decode(_dxSlider.getValue()), decode(_dtSlider.getValue()), lambdas, _wavepacket);
    }

	/***** MAIN FUNCTION *****/

    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    StolenCode frame = new StolenCode();
                    frame.setVisible(true);
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }
	
}
