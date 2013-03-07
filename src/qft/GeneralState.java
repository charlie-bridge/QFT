package qft;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import uk.ac.cam.cal56.maths.Complex;
import uk.ac.cam.cal56.qft.interactingtheory.Interaction;
import uk.ac.cam.cal56.qft.interactingtheory.State;

public class GeneralState implements State {

    // Class that represents a general state of the quantum system, and can
    // perform an integration step

    private int                      _systemSize; // Number of discrete spatial points in the system
    private double                   _epsilon;    // Distance between two adjacent points
    private double                   _mass;       // Mass of the boson (phi^2 interaction strength)
    private int                      _cutOffMom;  // Highest total Fock State momentum considered
    private int                      _numStates;  // Total number of Fock States considered (caution,
                                                   // index starts at 0)
    private int                      _phiPow;     // Power of the interaction term
    private InteractionMatrix        _interaction; // Interaction Matrix
    private double                   _dt;         // Time step size
    private double                   _time;       // Stores the system time
    private Complex[]                _coeffDerivs;       // Fock state coefficient derivatives at
                                                   // this point in time
    private Complex[]                _coeffs;     // Fock state coefficients at this point in time
    private double[]                 _frequencies;     // Frequencies of all of the states

    private Map<Interaction, Double> _lambdas;    // Interaction strengths

    private Complex[]                _nextCoeffs; // Required for integrator

    public GeneralState(int N, int Pmax, double m, double dx, double dt, Map<Interaction, Double> lambdas) {

        // Constructor for the general state

        // TODO Pmax not used by my general state but phi power and
        // cutoffMomentum needed, set to values for now

        _systemSize = N;
        _mass = m;
        _epsilon = dx;
        _cutOffMom = 50;
        _numStates = (((_cutOffMom * (_cutOffMom + 1)) / 2) + 1);
        _phiPow = 3;
        _dt = dt;

        _lambdas = new HashMap<Interaction, Double>();

        _interaction = new InteractionMatrix(_systemSize, _epsilon, _mass, _numStates, _phiPow);
        _coeffs = new Complex[_numStates];
        _coeffDerivs = new Complex[_numStates];
        _frequencies = new double[_numStates];

        // Calculate the interaction matrix

        _interaction.calcMatrix();

        // Calculate all the frequencies

        for (int i = 0; i < _numStates; i++) {

            _frequencies[i] = FockState.calcFrequency(i, _systemSize, _epsilon, _mass);

        }

        // Set the lambdas

        _lambdas = lambdas;

        // Set the state to the vacuum

        reset();

    }

    @Override
    public void setInteractionStrength(Interaction interaction, double lambda) {

        // Sets the interaction strengths

        _lambdas.put(interaction, lambda);

        // Perform the first Euler step again because the next coefficients are
        // no longer appropriate

        firstStep();

    }

    @Override
    public void step() {

        // Called in Carl's framework code, abstract away to allow multiple
        // integrators

        // step2ndOrderSymp();
        // simple();
        newIntegrator();

    }

    @Override
    public void step(int numSteps) {

        // Called in Carl's framework code, abstract away to allow multiple
        // integrators
        // Perform numSteps of steps

        for (int i = 0; i < numSteps; i++) {
            step();
        }

    }

    @SuppressWarnings("unused")
    private void step2ndOrderSymp() {

        // Second order symplectic integration step

        // TODO fix Complex addings

        Complex interactionPart = Complex.zero();

        // Increase the time

        _time += _dt;

        // Step the coefficient derivatives first

        for (int i = 0; i < _numStates; i++) {

            // First calculate the interaction part, multiply it by the
            // interaction strength

            for (Entry<Integer, Double> Hij : _interaction.getRow(i).entrySet()) {
                interactionPart = interactionPart.plus(_coeffs[Hij.getKey()].times(Hij.getValue()));

            }

            interactionPart.times(_lambdas.get(Interaction.PHI_CUBED));

            // Step the derivatives

            _coeffDerivs[i] = _coeffDerivs[i].plus(_coeffs[i].times(_frequencies[i]));
            _coeffDerivs[i] = _coeffDerivs[i].plus(interactionPart);
            _coeffDerivs[i] = _coeffDerivs[i].timesi(-1);
            _coeffDerivs[i] = _coeffDerivs[i].times(_lambdas.get(Interaction.PHI_CUBED));

            // Reset the interaction part

            interactionPart = Complex.zero();

        }

        // Now step the coefficients themselves

        for (int k = 0; k < _numStates; k++) {

            _coeffs[k] = _coeffs[k].plus(_coeffDerivs[k].times(_dt));

        }

    }

    @SuppressWarnings("unused")
    private void simple() {

        _time += _dt;

        Complex[] next = new Complex[_numStates];

        for (int i = 0; i < _numStates; i++) {

            Complex sum = Complex.zero();

            for (Entry<Integer, Double> Hij : _interaction.getRow(i).entrySet()) {
                sum = sum.plus(_coeffs[Hij.getKey()].times(Hij.getValue()));
            }

            sum = sum.times(_lambdas.get(Interaction.PHI_CUBED));
            sum = sum.plus(_coeffs[i].times(_frequencies[i]));
            Complex cdot = sum.timesi(-1);
            next[i] = _coeffs[i].plus(cdot.times(_dt));

        }

        _coeffs = next;

    }

    private void newIntegrator() {

        // Calculate the derivatives

        for (int i = 0; i < _numStates; i++) {

            _coeffDerivs[i] = Complex.zero();

            for (Entry<Integer, Double> Hij : _interaction.getRow(i).entrySet()) {
                _coeffDerivs[i] = _coeffDerivs[i].plus(_nextCoeffs[Hij.getKey()].times(Hij.getValue()));
            }

            _coeffDerivs[i] = _coeffDerivs[i].times(_lambdas.get(Interaction.PHI_CUBED));
            _coeffDerivs[i] = _coeffDerivs[i].plus(_nextCoeffs[i].times(_frequencies[i]));
            _coeffDerivs[i] = _coeffDerivs[i].timesi(-1);

        }

        // Store the _nextCoeffs for reassignment later

        Complex[] nextNextCoeffs = new Complex[_numStates];

        // Leapfrog

        for (int j = 0; j < _numStates; j++) {

            nextNextCoeffs[j] = _coeffs[j].plus(_coeffDerivs[j].times(2 * _dt));

        }

        // Set the old second to the new first

        _coeffs = _nextCoeffs;
        _nextCoeffs = nextNextCoeffs;

        // Step the time

        _time += _dt;

    }

    private void firstStep() {

        // Calculates a first Euler step

        Complex[] next = new Complex[_numStates];

        for (int i = 0; i < _numStates; i++) {

            _coeffDerivs[i] = Complex.zero();

            for (Entry<Integer, Double> Hij : _interaction.getRow(i).entrySet()) {
                _coeffDerivs[i] = _coeffDerivs[i].plus(_coeffs[Hij.getKey()].times(Hij.getValue()));
            }

            _coeffDerivs[i] = _coeffDerivs[i].times(_lambdas.get(Interaction.PHI_CUBED));

            _coeffDerivs[i] = _coeffDerivs[i].plus(_coeffs[i].times(_frequencies[i]));

            _coeffDerivs[i] = _coeffDerivs[i].timesi(-1);
            next[i] = _coeffs[i].plus(_coeffDerivs[i].times(_dt));

        }

        _nextCoeffs = next;

    }

    @Override
    public void reset(int... particleMomenta) {

        // Resets the general quantum state

        _time = 0.0;
        // _coeffs[0] = Complex.one();
        _coeffs[0] = Complex.zero();
        for (int i = 1; i < _numStates; i++) {
            _coeffs[i] = Complex.zero();
        }

        _coeffs[FockState.getIndex1PState(1, _systemSize, _epsilon, _mass)] = new Complex(0.70710678, 0);
        _coeffs[FockState.getIndex1PState(0, _systemSize, _epsilon, _mass)] = new Complex(0.5, 0);
        _coeffs[FockState.getIndex1PState(2, _systemSize, _epsilon, _mass)] = new Complex(0.5, 0);

        // Perform the first step

        firstStep();

    }

    @Override
    public double getTime() {

        // Returns the system time

        return _time;
    }

    @Override
    public Complex get0P() {

        // Gets the zero-particle state magnitude

        return _coeffs[0];
    }

    @Override
    public Complex[] get1PMom() {

        // Gets the one-particle state momenta magnitudes

        Complex[] onePMomMags = new Complex[_systemSize];
        Integer index = 0;

        // Set the stepping state to the first one momentum state

        for (int i = 0; i < _systemSize; i++) {

            // Find the index

            index = FockState.getIndex1PState(i, _systemSize, _epsilon, _mass);

            // Set the magnitude, if the state is over the cutoff set it to zero

            if ((index != null) && (index < _numStates)) {

                onePMomMags[i] = _coeffs[index];

            }
            else {

                onePMomMags[i] = Complex.zero();

            }

        }

        return onePMomMags;

    }

    @Override
    public Complex[][] get2PMom() {

        // Gets the two-particle momenta magnitudes

        Complex[][] twoPMomMags = new Complex[_systemSize][_systemSize];
        Integer index = 0;

        // Set the stepping states to the first two-particle state

        for (int i = 0; i < _systemSize; i++) {

            for (int j = i; j < _systemSize; j++) {

                // Calculate the index for this state

                index = FockState.getIndex2PState(i, j, _systemSize, _epsilon, _mass);

                // Store the magnitude, if the state is over the cutoff set it
                // to zero

                Complex value = Complex.zero();

                if ((index != null) && (index < _numStates)) {

                    value = _coeffs[index];
                    twoPMomMags[i][j] = value;

                }
                else {

                    twoPMomMags[i][j] = Complex.zero();

                }

                // If this does not lie on the diagonal, copy the value to its
                // mirror position too

                if (i != j) {
                    twoPMomMags[j][i] = value;
                }

            }

        }

        return twoPMomMags;

    }

    public double getModSquared() {

        // Returns the sum of the complex coefficient's moduli

        double returnVal = 0.0;

        for (int i = 0; i < _numStates; i++) {

            returnVal += _coeffs[i].modSquared();

        }

        return returnVal;
    }

    public double getRemainingProbability() {
        return 0.0;
    }

    public void setTimeStep(double dt) {

        _dt = dt;

    }

}
