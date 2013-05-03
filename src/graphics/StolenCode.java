package graphics;

import java.awt.Color;
import java.awt.EventQueue;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import javax.swing.JCheckBox;
import javax.swing.JSlider;

import qft.GeneralState;
import uk.ac.cam.cal56.graphics.Preset;
import uk.ac.cam.cal56.graphics.SimulatorFrame;
import uk.ac.cam.cal56.qft.Interaction;
import uk.ac.cam.cal56.qft.WavePacket;
import uk.ac.cam.cal56.qft.impl.MomentumWavePacket;

@SuppressWarnings("serial")
public class StolenCode extends SimulatorFrame {

	// Stolen QFTSandbox
	
    // @formatter:off
    protected String getFrameTitle()    { return "WOOOOO"; }
    protected Color getDisplayColor()   { return Color.BLACK; }
    protected Color getLabelColor()     { return Color.GRAY; }
    protected int getFrameWidth()       { return 1100; }
    protected int getFrameHeight()      { return 800; }
    protected int getPlotWidth()        { return 256; }
    protected int getPlotHeight()       { return getPlotWidth(); }
    protected int getNMin()             { return 2; }
    protected int getNMax()             { return 128; }
    protected int getPmaxMin()          { return 1; }
    protected int getPmaxMax()          { return 7; }
    protected double getDxMin()         { return 1.0e-3; }
    protected double getDxMax()         { return 10.0; }
    protected double getMMin()          { return 1.0e-3; }
    protected double getMMax()          { return 10.0; }
    protected double getDtMin()         { return 1.0e-5; }
    protected double getDtMax()         { return 1.0e-1; }
    protected int getStepsMin()         { return 1; }
    protected int getStepsMax()         { return 256; }
    protected double getLambdaMin()     { return 1.0e-7; }
    protected double getLambdaMax()     { return 1.0e2; }
    protected Preset getDefaultPreset() { return Preset.VACUUM; }
    // @formatter:on
	
	protected void setupQuantumState(WavePacket wavePacket) {
		Map<Interaction, Double> lambdas = new HashMap<Interaction, Double>();

		for (Entry<Interaction, JSlider> entry : interactionSliders.entrySet()) {
            Interaction interaction = entry.getKey();
            JSlider slider = entry.getValue();
            JCheckBox checkBox = interactionCheckBoxes.get(interaction);
            if (checkBox.isSelected())
                lambdas.put(interaction, decode(slider.getValue()));
		}

        int N = NSlider.getValue();
        if (wavePacket == null)
            wavePacket = new MomentumWavePacket(N);
        try {
            state = new GeneralState(N, PmaxSlider.getValue(), decode(mSlider.getValue()),
               decode(dxSlider.getValue()), decode(dtSlider.getValue()), lambdas, wavePacket);
            //state = new GeneralState(32, 5, 1e-1, 1e0, 1e-4, lambdas, wavePacket);
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
