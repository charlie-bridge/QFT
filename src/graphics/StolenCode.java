package graphics;

import java.awt.EventQueue;
import java.util.HashMap;
import java.util.Map;

import qft.GeneralState;
import uk.ac.cam.cal56.graphics.QFTSandbox;
import uk.ac.cam.cal56.qft.interactingtheory.Interaction;

@SuppressWarnings("serial")
public class StolenCode extends QFTSandbox {

	// Stolen QFTSandbox
	
	protected static final String FRAME_TITLE          = "BLEEERURUHJELAK";
	
	@Override
	protected void setupQuantumState() {
		Map<Interaction, Double> lambdas = new HashMap<Interaction, Double>();
		lambdas.put(Interaction.PHI_CUBED, decode(_lambdaSquaredSlider.getValue()));
        _state = new GeneralState(_NSlider.getValue(), _PmaxSlider.getValue(),
                decode(_mSlider.getValue()), decode(_dxSlider.getValue()), decode(_dtSlider.getValue()), lambdas);
        drawPlots();
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
