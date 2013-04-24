package graphics;

import java.awt.EventQueue;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import javax.swing.JCheckBox;
import javax.swing.JSlider;

import qft.GeneralState;
//import uk.ac.cam.cal56.graphics.SimulatorFrame;
import uk.ac.cam.cal56.qft.Interaction;
import uk.ac.cam.cal56.qft.WavePacket;
import uk.ac.cam.cal56.qft.impl.MomentumWavePacket;

@SuppressWarnings("serial")
public class StolenCode extends SimulatorFrameCh {

	// Stolen QFTSandbox
	
	protected static final String FRAME_TITLE = "BLEEERURUHJELAK";
	
	protected static final int                   FRAME_WIDTH           = 1240;
    protected static final int                   FRAME_HEIGHT          = 700;
	
	//private static final int _N_MAX = 256;
	
	protected void setupQuantumState(WavePacket wavePacket) {
		Map<Interaction, Double> lambdas = new HashMap<Interaction, Double>();

		for (Entry<Interaction, JSlider> entry : _interactionSliders.entrySet()) {
            Interaction interaction = entry.getKey();
            JSlider slider = entry.getValue();
            JCheckBox checkBox = _checkBoxes.get(interaction);
            if (checkBox.isSelected())
                lambdas.put(interaction, decode(slider.getValue()));
		}

        int N = _NSlider.getValue();
        if (wavePacket == null)
            wavePacket = new MomentumWavePacket(N);
        try {
            _state = new GeneralState(N, _PmaxSlider.getValue(), decode(_mSlider.getValue()),
               decode(_dxSlider.getValue()), decode(_dtSlider.getValue()), lambdas, wavePacket);
        }
        catch (IOException ioe) {
            System.out.println("Trouble above: " + ioe.getMessage());
        }
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
