package graphics;

import java.awt.EventQueue;

import qft.GeneralState;
import uk.ac.cam.cal56.graphics.QFTSandbox;

@SuppressWarnings("serial")
public class StolenCode extends QFTSandbox {

	// Stolen QFTSandbox
	
	@Override
	protected void setupQuantumState() {
        _state = new GeneralState(_N, _Pmax, _m, _dx, _dt, _lambda);
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
