package qft;

import uk.ac.cam.cal56.graphics.DensityPlot;
import uk.ac.cam.cal56.graphics.Plot;
import uk.ac.cam.cal56.maths.Complex;
import uk.ac.cam.cal56.maths.FFT;
import uk.ac.cam.cal56.maths.FourierTransform;
import uk.ac.cam.cal56.qft.interactingtheory.State;

public class GeneralState implements State {
	
	//Class that represents a general state of the quantum system, and can perform an integration step
	
	private int 				_systemSize;			//Number of discrete spatial points in the system
    private double 				_epsilon;				//Distance between two adjacent points
    private double 				_mass;					//Mass of the boson (phi^2 interaction strength)
    private int 				_cutOffMom;				//Highest total Fock State momentum considered
    private int 				_numStates;				//Total number of Fock States considered (caution, index starts at 0)
    private int 				_phiPow;				//Power of the interaction term
    private double 				_lambda;				//Interaction strength
    private InteractionMatrix 	_interaction;			//Interaction Matrix
    private double				_dt;					//Time step size
    private double 				_time;					//Stores the system time
    private Complex[]			_coeffDerivs;			//Fock state coefficient derivatives at this point in time
    private Complex[]           _coeffs;				//Fock state coefficients at this point in time
    private FourierTransform    _ft = new FFT();		//Fourier transform for calculating position probability magnitudes
    private double[] 			_frequencies;			//Frequencies of all of the states

	public GeneralState(int N, int Pmax, double m, double dx, double dt, double lambda){
		
		//Constructor for the general state
		
		//TODO Pmax not used by my general state but phi power and cutoffMomentum needed, set to values for now
		
		_systemSize = N;
		_mass = m;
		_epsilon = dx;
		_lambda = lambda;
		_cutOffMom = 200;
		_numStates = (((_cutOffMom*(_cutOffMom+1))/2)+1);
		_phiPow = 3;
		_dt = dt;
		
		_interaction = new InteractionMatrix(_systemSize, _epsilon, _mass, _numStates, _phiPow);
		_coeffs = new Complex[_numStates];
		_coeffDerivs = new Complex[_numStates];
		_frequencies = new double[_numStates];
		
		//Calculate the interaction matrix
		
		_interaction.calcMatrix();
		
		//Calculate all the frequencies
		
		for(int i=0; i<_numStates; i++) {
			
			_frequencies[i] = FockState.calcFrequency(i, _systemSize, _epsilon, _mass);
			
		}
		
		//Set the state to the vacuum
		
		reset();
		
	}
	
	@Override
	public void setInteractionStrength(double lambda) {
		
		//Sets the interaction strength lambda
		
		_lambda = lambda;
		
	}

	@Override
	public void step() {
		
		//Called in Carl's framework code, abstract away to allow multiple integrators
		
		step2ndOrderSymp();
		
	}
	
	private void step2ndOrderSymp() {
		
		//Second order symplectic integration step
		
		Complex interactionPart = Complex.zero();
		
		//Increase the time
		
		_time += _dt;
		
		//Step the coefficient derivatives first
		
		for(int i=0; i<_numStates; i++) {
			
			//First calculate the interaction part, multiply it by the interaction strength
			
			for(int j=0; j<_numStates; j++) {
				
				if (_interaction.getRow(j).get(i) != null) {
					interactionPart.plus(_coeffs[j].times(_interaction.getRow(j).get(i)));
				}

			}
			
			interactionPart.times(_lambda);
			
			//Step the derivatives
			
			_coeffDerivs[i].plus(_coeffs[i].times(_frequencies[i]));
			_coeffDerivs[i].plus(interactionPart);
			
			//Reset the interaction part
			
			interactionPart = Complex.zero();
			
		}
		
		//Now step the coefficients themselves
		
		for(int k=0; k<_numStates; k++) {
			
			_coeffs[k].plus(_coeffDerivs[k].times(_dt));
			
		}
		
	}

	@Override
	public void reset() {
		
		//Resets the general quantum state
		//RESETS TO VACUUM FOR NOW, RESET TO INITIAL LATER???
		
		_time = 0.0;
		_coeffs[0] = Complex.one();
		for(int i=1; i<_numStates; i++) {
			_coeffs[i] = Complex.zero();
		}
		
		//Set derivatives to zero too
		
		for(int j=1; j<_numStates; j++) {
			_coeffDerivs[j] = Complex.zero();
		}
		
		
	}

	@Override
	public double getTime() {
		
		//Returns the system time
		
		return _time;
	}

	@Override
	public double get0P() {
		
		//Gets the zero-particle state magnitude
		
		return _coeffs[0].modSquared();
	}

	@Override
	public double[] get1PMomenta() {
		
		//Gets the one-particle state momenta magnitudes
		
		double[] onePMomMags = new double[_systemSize];
		FockState stepState = new FockState(_systemSize, _epsilon, _mass);
		int stepIndex = 0;
		
		//Set the stepping state to the first one momentum state
		
		stepState.incrementState();

        for (int i=0; i<_systemSize; i++) {
        	
        	//Find the stepping state's index
        	
        	stepIndex = stepState.calcIndex();
        	
        	//Set the magnitude
        
        	onePMomMags[i] = _coeffs[stepIndex].modSquared();
        
        	//Increment the state
        	
        	stepState.incrementState();
        	
        }
        
        return onePMomMags;

	}

	@Override
	public double[] get1PPositions() {
		
		//Gets the one-particle state position magnitudes
		
		double[] onePPosMags = new double[_systemSize];
		Complex[] toTransform = new Complex[_systemSize];
		Complex[] transformed;
		FockState stepState = new FockState(_systemSize, _epsilon, _mass);
		int stepIndex = 0;
		
		//Prepare the list of coefficients to transform
		stepState.incrementState();
		
		for(int i=0; i<_systemSize; i++) {
			
			stepIndex = stepState.calcIndex();
			toTransform[i] = _coeffs[stepIndex];
			stepState.incrementState();
			
		}
		
		//Transform the list of coefficients
		
		transformed = _ft.transform(toTransform);
		
		//Set the magnitudes
		
		for (int j = 0; j<_systemSize; j++) {
			
			onePPosMags[j] = transformed[j].modSquared();
            
		}
		
		return onePPosMags;
	}

	@Override
	public double[][] get2PMomenta() {
		
		//Gets the two-particle momenta magnitudes
		
		double[][] twoPMomMags = new double[_systemSize][_systemSize];
		FockState stepState = new FockState(_systemSize, _epsilon, _mass);
		int stepIndex = 0;
		
		//Set the stepping states to the first two-particle state
		
		stepState.setCoeff(0, 2);
		
		for (int i = 0; i<_systemSize; i++) {
			
			for (int j = i; j<_systemSize; j++) {
				
				//Calculate the index for this state
				
				stepIndex = stepState.calcIndex();
				
				//Store the magnitude
				
				double value = _coeffs[stepIndex].modSquared();
				twoPMomMags[i][j] = value;
				
				//If this does not lie on the diagonal, copy the value to its mirror position too
				
				if (i != j) {
					twoPMomMags[j][i] = value;
				}
				
				//Increment the state
				
				stepState.incrementState();
				
			}
	            
		}
		
		return twoPMomMags;
	        
	}

	@Override
	public double[][] get2PPositions() {
		
		//Gets the two-particle position magnitudes
		
		double[][] twoPPosMags = new double[_systemSize][_systemSize];
		Complex[][] toTransform = new Complex[_systemSize][_systemSize];
		Complex[][] transformed;
		FockState stepState = new FockState(_systemSize, _epsilon, _mass);
		int stepIndex = 0;
		
		//Set the stepping state to the first two-particle state
		
		stepState.setCoeff(0, 2);

        for (int i=0; i<_systemSize; i++) {
        	
            for (int j=i; j <_systemSize; j++) {
            	
            	//Calculate the index for this state
            	
            	stepIndex = stepState.calcIndex();
            	
            	//Store the value that will be transformed
            	
                Complex value = _coeffs[stepIndex];
                toTransform[i][j] = value;
                
                //If this does not lie on the diagonal, copy the value to its mirror position too
                
                if (i != j) {
                    
                	toTransform[j][i] = value;
                	
                }
                
            }
            
        }
        
        //Transform the two-dimensional complex matrix
        
        transformed = _ft.transform2D(toTransform);
        
        for (int k = 0; k<_systemSize; k++) {

            for (int h = k; h<_systemSize; h++) {
            	
            	//Store the magnitude
            	
                double value = transformed[k][h].modSquared();
                twoPPosMags[k][h] = value;
                
                //If this does not lie on the diagonal, copy the value to its mirror position too
                
                if (k != h) {
                	
                	twoPPosMags[h][k] = value;
                	
                }
                
            }
            
        }
        
        return twoPPosMags;
        
	}

	@Override
	public void updatePlots(Plot p0m, Plot p0p, Plot p1m, Plot p1p,
			DensityPlot p2m, DensityPlot p2p) {
		// TODO Auto-generated method stub
		
	}
	

}
