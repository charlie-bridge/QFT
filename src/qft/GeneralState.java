package qft;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
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
    
    private InteractionMatrix        _interaction2; // Interaction Matrix for phi^2
    private InteractionMatrix        _interaction3; // Interaction Matrix for phi^3
    private InteractionMatrix        _interaction4; // Interaction Matrix for phi^4
    private double[]                 _frequencies;  // Frequencies of all of the states
    private Map<Interaction, Double> _lambdas;      // Interaction strengths
    
    private Complex[]                _coeffDerivs;  // Fock state coefficient derivatives at this point in time
    private Complex[]                _coeffs;       // Fock state coefficients at this point in time
    private Complex[]                _nextCoeffs;   // Required for leapfrog method
    
    private WavePacket               _wavePacket;   // User-set wavepacket
    
    private DataFile                 _dataFile;     // Data file to write to
    private String                   _fileName = "probConserv";  // Filename
    
    private double                   _1PSum;        // Used to calculate remaining probability
    private double                   _2PSum;        // Used to calculate remaining probability
    
    private List<Complex[]>          _eigenStates;  // Eigenstate storage
    private double                   _tau;          // Used in Power Iteration method

    public GeneralState(int N, int Pmax, double m, double dx, double dt, Map<Interaction, Double> lambdas, WavePacket wavePacket) throws IOException {

        // Constructor for the general state

        // TODO Pmax not used by my general state but cutoffMomentum needed, set to values for now

        _systemSize = N;
        _mass = m;
        _epsilon = dx;
        _cutOffMom = 40;
        _numStates = (((_cutOffMom * (_cutOffMom + 1)) / 2) + 1);
        _dt = dt;

        _interaction2 = new InteractionMatrix(_systemSize, _epsilon, _mass, _numStates, 2);
        _interaction3 = new InteractionMatrix(_systemSize, _epsilon, _mass, _numStates, 3);
        _interaction4 = new InteractionMatrix(_systemSize, _epsilon, _mass, _numStates, 4);
        _coeffs = new Complex[_numStates];
        _nextCoeffs = new Complex[_numStates];
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

        // Set the state

        setWavePacket(wavePacket);
        
        //Set up file to write to
        
        _dataFile = new DataFile(_fileName);  
        
        //Set up eigenstate storage
        
        _eigenStates = new ArrayList<Complex[]>();
        
        
        //testPhi2AsMass();

    }

    @Override
    public void setInteractionStrength(Interaction interaction, double lambda) {

        // Sets the interaction strengths

        if(_lambdas.containsKey(interaction)) {
            _lambdas.put(interaction, lambda);
        }

        // Perform the first Euler step again because the next coefficients are no longer valid

        firstStep();

    }

    @Override
    public void step() {

        // Called in Carl's framework code, abstract away to allow multiple
        // integrators

        //semiImplicitEuler();
        leapfrog();
        //nystrom();
        //firstEuler();
        
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

    private Complex[] applyHamiltonian(Complex[] vector, boolean withMinusI) {
        
        //Applies the Hamiltonian to the passed vector of coefficients
        
        Complex[] returnVector = new Complex[_numStates];
        
        for (int i = 0; i < _numStates; i++) {

            returnVector[i] = Complex.zero();

            if(_lambdas.containsKey(Interaction.PHI_SQUARED)) {
                for (Entry<Integer, Double> Hij : _interaction2.getRow(i).entrySet()) {
                    
                    returnVector[i] = returnVector[i].plus(vector[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_SQUARED)));
                    
                }
            }
            
            if(_lambdas.containsKey(Interaction.PHI_CUBED)){
                for (Entry<Integer, Double> Hij : _interaction3.getRow(i).entrySet()) {
                    
                    returnVector[i] = returnVector[i].plus(vector[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_CUBED)));
                    
                }
            }
            
            if(_lambdas.containsKey(Interaction.PHI_FOURTH)) {
                for (Entry<Integer, Double> Hij : _interaction4.getRow(i).entrySet()) {
                        
                    returnVector[i] = returnVector[i].plus(vector[Hij.getKey()].times(Hij.getValue()*_lambdas.get(Interaction.PHI_FOURTH)));
                    
                }
            }
            
            returnVector[i] = returnVector[i].plus(vector[i].times(_frequencies[i])); 
            
            //Check for multiplying being minus i
            
            if(withMinusI == true) {
                returnVector[i] = returnVector[i].timesi(-1);
            }
            
        }
        
        return returnVector;
        
    }

    private void semiImplicitEuler() {
        
        //Cpply the semi-implicit Euler method integration step
        
        Complex[] firstHalf = new Complex[_numStates];
        Complex[] secondHalf = new Complex[_numStates];
        
        firstHalf = applyHamiltonian(_coeffs, true);
        secondHalf = applyHamiltonian(firstHalf, true);
        
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
        
        addor = applyHamiltonian(_nextCoeffs, true);
        
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
        
        firstHalf = applyHamiltonian(_nextCoeffs, true);
        secondHalf = applyHamiltonian(_coeffs, true);
        
        for(int i=0; i<_numStates; i++) {
            
            _nextCoeffs[i] = _coeffs[i].plus(firstHalf[i].times(0.5*_dt));
            _nextCoeffs[i] = _nextCoeffs[i].plus(secondHalf[i].times(1.5*_dt));
            
        }
        
        _coeffs = store;
        
        // Step the time

        _time += _dt;
    }
    
    private void firstEuler() {
        
        _coeffDerivs = applyHamiltonian(_coeffs, true);
        
        for(int i=0; i<_numStates; i++) {
            
            _coeffs[i] = _coeffs[i].plus(_coeffDerivs[i].times(_dt));
            
        }
        
        _time += _dt;
        
    }
    
    private void newIntegrator() {

        //OBSOLETE
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
        
        //Perform a first order Euler step to calculate the second set of coefficients for the leapfrog method
        
        _coeffDerivs = applyHamiltonian(_coeffs, true);
        
        for(int i=0; i<_numStates; i++) {
            
            _nextCoeffs[i] = _coeffs[i].plus(_coeffDerivs[i].times(_dt));
            
        }
        
        _time += _dt;
    }
    
    private void firstStepOLD() {

        //OBSOLETE
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

        // Gets the zero-particle state coefficient

        return _coeffs[0];
                
    }

    @Override
    public Complex[] get1PMom() {

        _1PSum = 0.0;
        
        // Gets the one-particle state momenta coefficients

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

        // Gets the two-particle momenta coefficients

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

        // Returns the sum of the complex coefficient's moduli squared

        double returnVal = 0.0;
        
        FockState state = new FockState(_systemSize, _mass, _epsilon);

        for (int i = 0; i < _numStates; i++) {
               
            state.setAsIndex(i);
            
            returnVal += (_coeffs[i].modSquared());

        }

        return returnVal;
    }

    public Double getRemainingProbability() {
        
        double sum012;
        sum012 = _coeffs[0].modSquared() + _1PSum + _2PSum;
        
        return (1-sum012);
        
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
    
    public double getTotalEnergy(){
        
        return getTotalStateEnergy(null);
        
    }
    
    private double getTotalStateEnergy(Complex[] state) {
        
        //Returns the state's total energy, or a specific state if one is supplied
        
        double energy = 0.0;
        double imagPart = 0.0;
        
        Complex[] hammedVector = new Complex[_numStates];
        
        if(state == null) {
            hammedVector = applyHamiltonian(_coeffs, false);
        }
        else {
            hammedVector = applyHamiltonian(state, false);
        }
        
        for(int i=0; i<_numStates; i++) {
            
            if(state == null) {
                energy += ((_coeffs[i].conj()).times(hammedVector[i])).real();
            }
            else {
                energy += ((state[i].conj()).times(hammedVector[i])).real();
            }
        }
        
        //TODO
        //double zp = 0.0;
        //for(int i=0; i<_systemSize; i++) {
        //    zp += FockState.calcFrequency(i, _systemSize, _epsilon, _mass);
        //}
        //zp=(zp/2);
        
        //return (energy +zp);
        
        return energy;
        
    }
    
    private void writeToFile() throws IOException {
        
        //Writes something to file
        
        Double[] line = new Double[3];
        line[0] = _time;
        line[1] = getModSquared();
        line[2] = getTotalStateEnergy(null);
        
        _dataFile.writeLine(line);
        
    }
    
    private Complex[] findState(int stateNumber) {
        
        //Uses the Power Law method to find an eigenstate state of the full Hamiltonian
        //For stateNumber=0 find the ground state
        //For stateNumber>0 use the previously found states to find the stateNumber th eigenstate
        
        //Find the largest eigenvaue to ensure fastest convergence
        
        if(stateNumber == 0) {
            double largest = findLargestEval();
            _tau = Math.ceil(largest);
        }
        
        //Set the initial vector
        
        Complex[] vector1 = new Complex[_numStates];
        Complex[] vector2 = new Complex[_numStates];
        for(int i=0; i<_numStates; i++) {
            if(i<10) {
                vector1[i] = Complex.one().divide(Math.sqrt(10.0));
            }
            else{
                vector1[i] = Complex.zero();
            }
        }
        
        //Set the tolerance value and maximum number of iterations
        
        double tolerance = 1e-11;
        int maxIterations = 100000;
        double norm = 0.0;
        int number = 0;
        double difference;
        double eigenVecAmount;
        boolean flag = true;
        
        //Repeat the algorithm until the tolerance is met or the maximum number of iterations is hit
        
        while(flag && (number<maxIterations)) {
              
            //Remove all the eigenvector components of the state corresponding to
            //eigenvectors found with eigenvalues lower than this one
            
            for(int j=0; j<stateNumber; j++){
                
                eigenVecAmount = dotProduct(_eigenStates.get(j), vector1);
                
                for(int k=0; k<_numStates; k++) {
                    
                    vector1[k] = vector1[k].minus((_eigenStates.get(j))[k].times(eigenVecAmount)); 
                    
                }
                
                //Renormalise each time
                
                norm = normVector(vector1, null);
                
                for(int k=0; k<_numStates; k++) {
                    
                    vector1[k] = vector1[k].divide(norm);
                    
                }
                
            }      
            
            //Apply the algorithm
            
            vector2 = applyHamiltonian(vector1, false);
            
            for(int i=0; i<_numStates; i++) {
                vector2[i] = (vector1[i].times(_tau)).minus(vector2[i]);
            }
            
            norm = normVector(vector2, null);
            
            for(int i=0; i<_numStates; i++) {
                vector2[i] = vector2[i].divide(norm);
            }
            
            //Determine the change in vector caused by the algorithm for tolerance check
            
            difference = normVector(vector1, vector2);
            
            if(difference<tolerance) {
                flag = false;
            }
            
            //Prepare for next iteration
            
            vector1 = vector2;
            
            number++;
            
            
            if(number%1000 == 0) {
                System.out.println(number);
                System.out.println(vector2[0] + " " + vector2[1] + " " + vector2[2] + " " + vector2[3] + " " + vector2[4] + " " + vector2[5]);
            }
            
        }
        
        //Add this eigenstate to the stored list and return it
        
        _eigenStates.add(stateNumber, vector2);        
        
        
        System.out.println(vector2[0] + " " + vector2[1] + " " + vector2[2] + " " + vector2[3] + " " + vector2[4] + " " + vector2[5]);
        
        
        return vector2;
        
    }
    
    private double findLargestEval() {
        
        //Finds the largest eigenvalue of the full Hamiltonian using the simple power iteration method
        
        Complex[] vector1 = new Complex[_numStates];
        Complex[] vector2 = new Complex[_numStates];
        
        for(int i=0; i<_numStates; i++) {
           vector1[i] = Complex.zero();
        }
        
        vector1[_numStates-1] = Complex.one();
        
        double tolerance = 0.00000000001;
        double norm = 0.0;
        int number = 0;
        double difference;
        boolean flag = true;
        
        while(flag && (number<100000)) {
            
            vector2 = applyHamiltonian(vector1, false);
            
            norm = normVector(vector2, null);
            
            for(int i=0; i<_numStates; i++) {
                vector2[i] = vector2[i].divide(norm);
            }
            
            difference = normVector(vector1, vector2);
            
            if(difference<tolerance) {
                flag = false;
            }
            
            vector1 = vector2;
            
            number++;
            
        }
        
        return normVector(applyHamiltonian(vector2, false), null);
        
    }
    
    private double dotProduct(Complex[] vec1, Complex[] vec2) {
        
        //Returns the dot product of two vectors
        
        double returnval = 0.0;
        Complex returnComp = Complex.zero();
        
        for(int i=0; i<_numStates; i++) {
            returnComp = returnComp.plus(((vec1[i].conj()).times(vec2[i])));
        }
        
        //Make sure the sign is not lost
        
        int sign = (int)(returnComp.real()/Math.abs(returnComp.real()));
        returnval = ((returnComp.mod())*(double)(sign));
        return returnval;
        
    }
    
    public void setToGroundState() {
        
        //Sets the system to the ground state
        
        _coeffs = findState(0);
        
    }
    
    private double normVector(Complex[] vector1, Complex[] vector2) {
        
        //Returns the norm of the vector1 if no vector2 is supplied
        //Returns the norm of the difference between the vectors if both are supplied
        
        double returnValue = 0.0;
        
        if(vector2 == null) {
            for(int i=0; i<_numStates; i++) {
                returnValue += vector1[i].modSquared();
            }
        }
        else {
            for(int i=0; i<_numStates; i++) {
                returnValue += (vector1[i].minus(vector2[i])).modSquared();
            }
        }
        
        returnValue = Math.sqrt(returnValue);
        
        return returnValue;
        
    }
    
    private void writeBasisEnergies() throws IOException {
        
        //Writes the basis state energies to file in order of index
        
        Complex[] basisState = new Complex[_numStates];
        Double[] energy = new Double[1];
        
        for(int i=0; i<_numStates; i++) {
            
            for(int j=0; j<_numStates; j++) {
                basisState[j] = Complex.zero();            
            }
            basisState[i] = Complex.one();
            
            energy[0] = getTotalStateEnergy(basisState);
            
            _dataFile.writeLine(energy);
       
        }
        
    }

    private void testPhi2AsMass() throws IOException {
        
        //Produces the ground state energy as a function of lambda2
        
        //Set lambda3 and lambda4 to zero
        
        _lambdas.put(Interaction.PHI_CUBED, null);
        _lambdas.put(Interaction.PHI_FOURTH, null);
        _interaction2.calcMatrix();
        
        //For each value of _lambda compute the ground state energy and print
        
        Double data[] = new Double[2];
        
        for(int i=0; i<7; i++) {
            
            //data[0] = Math.pow(10.0, (i-5));
            data[0] = Math.pow(10.0, -2);
            _lambdas.put(Interaction.PHI_SQUARED, data[0]);
            
            for(int x=0; x<10; x++) {
                _coeffs = findState(x);   
                data[1] = getTotalStateEnergy(_coeffs);
                _dataFile.writeLine(data);
            }
            
        }
        

    }

}
