package qft;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import uk.ac.cam.cal56.maths.Complex;
import uk.ac.cam.cal56.qft.Interaction;
import uk.ac.cam.cal56.qft.State;
import uk.ac.cam.cal56.qft.WavePacket;
import utilities.DataFile;

public class GeneralState implements State {

    // Class that represents a general state of the quantum system, and can
    // perform an integration step

    private int                      _systemSize;   // Number of discrete spatial points in the system
    private double                   _epsilon;      // Distance between two adjacent points
    private double                   _mass;         // Mass of the boson (phi^2 interaction strength)
    private int                      _cutOffMom;    // Highest total Fock State momentum considered
    private int                      _numStates;    // Total number of Fock States considered (caution, index starts at 0)
    private double                   _dt;           // Time step size
    private double                   _time;         // Stores the system time
    //private int                      _phiPow;       // Power of the interaction term
    
    private InteractionMatrix        _interaction2; // Interaction Matrix for phi^2
    private InteractionMatrix        _interaction3; // Interaction Matrix for phi^3
    private InteractionMatrix        _interaction4; // Interaction Matrix for phi^4
    private double[]                 _frequencies;  // Frequencies of all of the state
    private Map<Interaction, Double> _lambdas;      // Interaction strengths
    
    private Complex[]                _coeffDerivs;  // Fock state coefficient derivatives at this point in time
    private Complex[]                _coeffs;       // Fock state coefficients at this point in time
    private Complex[]                _nextCoeffs;   // Required for integrator
    
    private WavePacket               _wavePacket;   // User-set wavepacket
    
    private DataFile                 _dataFile;     // Data file to write to
    
    private double                   _1PSum;        // Used to calculate remaining probability
    private double                   _2PSum;        // Used to calculate remaining probability

    public GeneralState(int N, int Pmax, double m, double dx, double dt, Map<Interaction, Double> lambdas, WavePacket wavePacket) throws IOException {

        // Constructor for the general state

        // TODO Pmax not used by my general state but cutoffMomentum needed, set to values for now

        _systemSize = N;
        _mass = m;
        _epsilon = dx;
        _cutOffMom = 40;
        _numStates = (((_cutOffMom * (_cutOffMom + 1)) / 2) + 1);
        _dt = dt;
        
        //_phiPow = 3;

        _lambdas = new HashMap<Interaction, Double>();

        _interaction2 = new InteractionMatrix(_systemSize, _epsilon, _mass, _numStates, 2);
        _interaction3 = new InteractionMatrix(_systemSize, _epsilon, _mass, _numStates, 3);
        _interaction4 = new InteractionMatrix(_systemSize, _epsilon, _mass, _numStates, 4);
        _coeffs = new Complex[_numStates];
        _coeffDerivs = new Complex[_numStates];
        _frequencies = new double[_numStates];
        
        // Set the lambdas

        _lambdas = lambdas;

        // Calculate the interaction matrix(ces), but only if necessary

        if(_lambdas.get(Interaction.PHI_SQUARED) != null) {
            _interaction2.calcMatrix();
        }
        if(_lambdas.get(Interaction.PHI_CUBED) != null) {
            _interaction3.calcMatrix();
        }
        if(_lambdas.get(Interaction.PHI_FOURTH) != null) {
            _interaction4.calcMatrix();
        }

        // Calculate all the frequencies

        FockState energyState = new FockState(N, dx, m);
        for (int i = 0; i < _numStates; i++) {

            energyState.setAsIndex(i);
            _frequencies[i] = energyState.calcEnergy();

        }

        // Set the state to the vacuum

        setWavePacket(wavePacket);
        
        //Set up file to write to
        
        _dataFile = new DataFile("test");
        

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

        //newIntegrator();
        //semiImplicitEuler();
        //leapfrog();
        //nystrom();
        firstEuler();
        
        try {
            if((((int)(_time/_dt))%1000) == 0) {
                writeToFile();
            }
        } catch (IOException ioe) {
            System.out.println("Trouble: " + ioe.getMessage());
        }

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

    private Complex[] applyHamiltonian(Complex[] vector) {
        
        //Applies the Hamiltonian to the passed vector of coefficients
        
        Complex[] returnVector = new Complex[_numStates];
        
        for (int i = 0; i < _numStates; i++) {

            returnVector[i] = Complex.zero();

            if(_lambdas.get(Interaction.PHI_SQUARED) != null) {
                for (Entry<Integer, Double> Hij : _interaction2.getRow(i).entrySet()) {
                    
                    returnVector[i] = returnVector[i].plus(vector[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_SQUARED)));
                    
                }
            }
            
            if(_lambdas.get(Interaction.PHI_CUBED) != null) {
                for (Entry<Integer, Double> Hij : _interaction3.getRow(i).entrySet()) {
                    
                    returnVector[i] = returnVector[i].plus(vector[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_CUBED)));
                    
                }
            }
            
            if(_lambdas.get(Interaction.PHI_FOURTH) != null) {
                for (Entry<Integer, Double> Hij : _interaction4.getRow(i).entrySet()) {
                    
                    returnVector[i] = returnVector[i].plus(_nextCoeffs[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_FOURTH)));
                    
                }
            }
            
            returnVector[i] = returnVector[i].plus(vector[i].times(_frequencies[i])); 
            
            //MAYBE DONT HAVE THE PHASE HERE
            
            returnVector[i] = returnVector[i].timesi(-1);
            
        }
        
        return returnVector;
        
    }

    private void semiImplicitEuler() {
        
        //Cpply the semi-implicit Euler method integration step
        
        Complex[] firstHalf = new Complex[_numStates];
        Complex[] secondHalf = new Complex[_numStates];
        
        firstHalf = applyHamiltonian(_coeffs);
        secondHalf = applyHamiltonian(firstHalf);
        
        for(int i=0; i<_numStates; i++) {
            
            _coeffs[i] = _coeffs[i].plus(firstHalf[i].times(_dt));
            _coeffs[i] = _coeffs[i].plus(secondHalf[i].times(_dt*_dt));
            
        }
        
        // Step the time

        _time += _dt;
        
    }
    
    private void leapfrog() {
        
        //Apply the leapfrog method integration step
        
        Complex[] addor = new Complex[_numStates];
        Complex[] store = new Complex[_numStates];
        
        for(int i=0; i<_numStates; i++) {
            store[i] = _nextCoeffs[i].times(1.0);
        }
        
        addor = applyHamiltonian(_nextCoeffs);
        
        for(int i=0; i<_numStates; i++) {
            
            _nextCoeffs[i] = _coeffs[i].plus(addor[i].times(2*_dt));
            
        }
        
        _coeffs = store;
        
        // Step the time

        _time += _dt;
        
    }
    
    private void nystrom() {
        
        Complex[] firstHalf = new Complex[_numStates];
        Complex[] secondHalf = new Complex[_numStates];
        Complex[] store = new Complex[_numStates];
        
        for(int i=0; i<_numStates; i++) {
            store[i] = _nextCoeffs[i].times(1.0);
        }
        
        firstHalf = applyHamiltonian(_nextCoeffs);
        secondHalf = applyHamiltonian(_coeffs);
        
        for(int i=0; i<_numStates; i++) {
            
            _nextCoeffs[i] = _coeffs[i].plus(firstHalf[i].times(0.5*_dt));
            _nextCoeffs[i] = _nextCoeffs[i].plus(secondHalf[i].times(1.5*_dt));
            
        }
        
        _coeffs = store;
        
        // Step the time

        _time += _dt;
    }
    
    private void firstEuler() {
        
        _coeffDerivs = applyHamiltonian(_coeffs);
        
        for(int i=0; i<_numStates; i++) {
            
            _coeffs[i] = _coeffs[i].plus(_coeffDerivs[i].times(_dt));
            
        }
        
        _time += _dt;
        
    }
    
    private void newIntegrator() {

        // Calculate the derivatives

        for (int i = 0; i < _numStates; i++) {

            _coeffDerivs[i] = Complex.zero();

            if(_lambdas.get(Interaction.PHI_SQUARED) != null) {
                for (Entry<Integer, Double> Hij : _interaction2.getRow(i).entrySet()) {
                    
                    _coeffDerivs[i] = _coeffDerivs[i].plus(_nextCoeffs[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_SQUARED)));
                    
                }
            }
            
            if(_lambdas.get(Interaction.PHI_CUBED) != null) {
                for (Entry<Integer, Double> Hij : _interaction3.getRow(i).entrySet()) {
                    
                    _coeffDerivs[i] = _coeffDerivs[i].plus(_nextCoeffs[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_CUBED)));
                    
                }
            }
            
            if(_lambdas.get(Interaction.PHI_FOURTH) != null) {
                for (Entry<Integer, Double> Hij : _interaction4.getRow(i).entrySet()) {
                    
                    _coeffDerivs[i] = _coeffDerivs[i].plus(_nextCoeffs[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_FOURTH)));
                    
                }
            }
            
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

            if(_lambdas.get(Interaction.PHI_SQUARED) != null) {
                for (Entry<Integer, Double> Hij : _interaction2.getRow(i).entrySet()) {
                    
                    _coeffDerivs[i] = _coeffDerivs[i].plus(_coeffs[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_SQUARED)));
                    
                }
            }
            
            if(_lambdas.get(Interaction.PHI_CUBED) != null) {
                for (Entry<Integer, Double> Hij : _interaction3.getRow(i).entrySet()) {
                    
                    _coeffDerivs[i] = _coeffDerivs[i].plus(_coeffs[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_CUBED)));
                    
                }
            }
            
            if(_lambdas.get(Interaction.PHI_FOURTH) != null) {
                for (Entry<Integer, Double> Hij : _interaction3.getRow(i).entrySet()) {
                    
                    _coeffDerivs[i] = _coeffDerivs[i].plus(_coeffs[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_FOURTH)));
                    
                }
            }
            
            _coeffDerivs[i] = _coeffDerivs[i].plus(_coeffs[i].times(_frequencies[i]));

            _coeffDerivs[i] = _coeffDerivs[i].timesi(-1);
            
            next[i] = _coeffs[i].plus(_coeffDerivs[i].times(_dt));

        }

        _nextCoeffs = next;

    }

    @Override
    public double getTime() {

        // Returns the system time

        return _time;
    }

    @Override
    public Complex getVacuum() {

        // Gets the zero-particle state magnitude

        return _coeffs[0];
                
    }

    @Override
    public Complex[] get1PMom() {

        _1PSum = 0.0;
        
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
                _1PSum += _coeffs[index].modSquared();
         
            }
            else {

                onePMomMags[i] = Complex.zero();

            }

        }

        return onePMomMags;

    }

    @Override
    public Complex[][] get2PMom() {
        
        _2PSum = 0.0;

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
                    _2PSum += _coeffs[index].modSquared();

                }
                else {

                    twoPMomMags[i][j] = Complex.zero();

                }

                // If this does not lie on the diagonal, copy the value to its
                // mirror position too
                //TODO DIVIDE BY 2????????????????

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
        
        FockState state = new FockState(_systemSize, _mass, _epsilon);

        for (int i = 0; i < _numStates; i++) {
            
            
            state.setAsIndex(i);
            
            
            //returnVal += (_coeffs[i].modSquared() * (Math.pow(state.norm(), 2.0)) * 2 * Math.PI * (_epsilon/_systemSize)));;
            
            
            returnVal += (_coeffs[i].modSquared());

        }

        return returnVal;
    }

    public Double getRemainingProbability() {
        
        double sum012;
        sum012 = _coeffs[0].modSquared() + _1PSum + _2PSum;
        
        return (1-sum012)
                ;
    }

    public void setTimeStep(double dt) {

        _dt = dt;

    }

    @Override
    public void reset() {
        
        // Reset stuff
        
        _coeffs = CoeffAdaptor.setCoeffs(_systemSize, _numStates, _wavePacket.getCoefficients((_systemSize + 1) * (_systemSize + 2) / 2));
        
        _time = 0.0;
        firstStep();
        
    }

    @Override
    public void setWavePacket(WavePacket wavePacket) {
        
        // Reset stuff
        
        _wavePacket = wavePacket;
        reset();
        
    }

    @Override
    public int getN() {
        
        // Gets the system size
        
        return _systemSize;
    }
    
    private double getTotalEnergy() {
        
        //Returns the state's total energy
        
        double energy = 0.0;
        
        for(int i=0; i<_numStates; i++) {
            
            energy += (_coeffs[i].modSquared()*_frequencies[i]);
            
        }
        
        return energy;
        
    }
    
    private void writeToFile() throws IOException {
        
        //Writes something to file
        
        Double[] line = new Double[2];
        line[0] = _time;
        line[1] = getTotalEnergy();
        _dataFile.writeLine(line);
        
    }
    
}
