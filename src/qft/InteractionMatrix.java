package qft;

import java.util.*;

import utilities.Factorial;

public class InteractionMatrix {
    
    //Class for calculating and storing an interaction matrix for theory of general interaction power
    
    private int _systemSize;							//Number of discrete spatial points in the system
    private double _epsilon;							//Distance between two adjacent points
    private double _mass;								//Mass of the boson (phi^2 interaction strength)
    private int _cutOffMom;								//Highest total Fock State momentum considered
    private int _numStates;								//Total number of Fock States considered
    private int _phiPow;								//Power of the interaction term
    private int _rightMostP;							//Farthest right operator momentum
    private double _operatorFactors;					//Factors resulting from application of operators
    private List<Map<Integer, Double>> _elements;		//Store of all the matrix elements
    
    public InteractionMatrix(int systemSize, double epsilon, double mass, int cutOffMom, int phiPow) {
        
        //Cut off momentum should be greater than or equal to one
    	
    	//TODO Test declaring Map and Lists sizes against not
        
        _systemSize = systemSize;
        _epsilon = epsilon;
        _mass = mass;
        _cutOffMom = cutOffMom;
        _phiPow = phiPow;
        _numStates = (((_cutOffMom*(_cutOffMom+1))/2) + 1);
        _elements = new ArrayList<Map<Integer, Double>>(_numStates);
        
    }
    
    public void calcMatrix() {
        
        FockState rowState = new FockState(_systemSize, _epsilon, _mass);
        FockState opedState = new FockState(_systemSize, _epsilon, _mass);
        int column;
        
        int[] momenta = new int[_phiPow-1];
        
        for(int row=0; row<_numStates; row++) {
            
            //Run through the matrix rows
        	
        	for(int i=0; i<(_phiPow-1); i++) {
                momenta[i]=0;
            }
            
            rowState.setAsIndex(row);
            
            for(int j=0; j<((int)Math.pow(2, (_phiPow - 1))); j++) {
                
                //Run through the various momenta combinations
            
                for(int k=0; k<((int)Math.pow(2, _phiPow)); k++) {
                    
                    //Run through the various operator combinations
                    
                    opedState.makeSameAs(rowState);
                    applyOperators(opedState, momenta, k);
                    
                    if(opedState.isValid()==true) {
                        
                        //If the resulting state is valid, store the result
                        
                        column = opedState.calcIndex();
                        store(row, column, momenta);
                        
                    }
                
                }
            
                //Increment the momenta label array
                
                incrementMomentaLabel(momenta);
                
            }
            
        }
        
    }
    
    public void incrementMomentaLabel(int[] momenta) {
    	
    	//PUBLIC FOR TESTING ONLY, DO NOT USE OUTSIDE THIS CLASS
        
        //Increments the momenta label array
        //Add one to first value, but if that is = (_systemSize-1) then set it to zero
        //increase the next one, but if that is = (_systemSize-1) etc.
        
        boolean keepGoing = true;
        
        for(int i=0; (i<(_phiPow-1) && (keepGoing == true)); i++) {
            
            if(momenta[i]<(_systemSize-1)) {
                momenta[i]++;
                keepGoing = false;
            }
            else{
                momenta[i] = 0;
            }
            
        }
        
        //If keepGoing is still true, all values must be N, so all combinations have been considered
        //WATCH OUT FOR THIS
        
    }
    
    public void applyOperators(FockState passedRow, int[] momenta, int opTypes) {
    	
    	//PUBLIC FOR TESTING ONLY, DO NOT USE OUTSIDE THIS CLASS
        
        //Applies operators to the passed row FockState according to the integer i
        //We have 2^phiPow possibilities for the operators - so use binary
        //For each 1 bit use annihilation, for each 0 bit use creation (on the passed <bra|)
        //The rightmost operator is special, its momentum is set by the others
    	//this is stored outside this method for use later
        
        _rightMostP = 0;
        _operatorFactors = 1.0;
        
        for(int i=0; i<(_phiPow-1); i++) {
        	
        	//Run through the operators from left to right
            
            if(((opTypes & (1 << i)) != 0)) {
            	
            	//This bit is a 1, apply an annihilation operator to the relevant momentum mode
            	//Update the running operator factor and rightmost operator momentum
            	//If this operation invalidates the state set the factor to zero
                
                if(passedRow.getCoeff(momenta[i]) > 0) {
            		_operatorFactors = (_operatorFactors * Math.sqrt(passedRow.getCoeff(momenta[i])));
            	}
                else {
                	_operatorFactors = 0.0;
                }
                
                passedRow.applyAnnihil(momenta[i]);
                _rightMostP = _rightMostP - momenta[i];
                
            }
            else {
            	
            	//This bit is a 0, apply a creation operator to the relevant momentum mode
            	//Update the running operator factor and rightmost operator momentum
                
            	_operatorFactors = (_operatorFactors * Math.sqrt(passedRow.getCoeff(momenta[i]) + 1.0));
                passedRow.applyCreation(momenta[i]);
                _rightMostP = _rightMostP + momenta[i];
                
            }
            
        }
        
        //Put the rightmost momenta in the required range before using it
        
        if(((opTypes & (1 << (_phiPow-1))) != 0)) {
        	
        	//This bit is a 1, apply an annihilation operator to the rightmost operator momentum mode
        	//Complete the operator factor for this set of operators
        	//If this operation invalidates the state set the factor to zero
        	
        	if (_rightMostP >= 0) {
        		_rightMostP = (Math.abs(_rightMostP % _systemSize));
        	}
        	else {
        		_rightMostP = _systemSize - (Math.abs(_rightMostP % _systemSize));
        	}
            
        	if(passedRow.getCoeff(_rightMostP) > 0) {
        		_operatorFactors = (_operatorFactors * Math.sqrt(passedRow.getCoeff(_rightMostP)));
        	}
        	else {
        		_operatorFactors = 0.0;
        	}
        	
        	passedRow.applyAnnihil(_rightMostP);
            
        }
        else {
        	
        	//This bit is a 0, apply a creation operator to the rightmost operator momentum mode
        	//Complete the operator factor for this set of operators
        	
        	if (_rightMostP >= 0) {
        		_rightMostP = _systemSize - (Math.abs(_rightMostP % _systemSize));
        	}
        	else {
        		_rightMostP = (Math.abs(_rightMostP % _systemSize));
        	}
        	
        	_rightMostP = Math.abs(_rightMostP % _systemSize);
        	_operatorFactors = (_operatorFactors * Math.sqrt(passedRow.getCoeff(_rightMostP) + 1.0));
            passedRow.applyCreation(_rightMostP);
            
            
        }
        
    }
    
    public void store(int row, int column, int[] momenta) {
    	
    	//PUBLIC FOR TESTING ONLY, DO NOT USE OUTSIDE THIS CLASS
        
    	//Calculate the value to store
    	//The rightmost momentum value has already been calculated by applyOperators
    	
    	double maElVal = (_epsilon * Math.pow((_systemSize*_epsilon * 2.0), -(_phiPow/2.0)) * (1/Factorial.calc(_phiPow)));
    	maElVal = (maElVal * _operatorFactors);
    	maElVal = (maElVal * Math.pow(FockState.calcFrequency(_rightMostP, _systemSize, _epsilon, _mass), -0.5));
    	
    	for(int i=0; i<(_phiPow-1); i++) {
    		
    		maElVal = (maElVal * Math.pow(FockState.calcFrequency(momenta[i], _systemSize, _epsilon, _mass), -0.5));
    		
    	}
    	
    	//Check if a value already exists at this location, if it does, just add to it
    	
    	if(_elements.get(row).containsKey(column)) {
    		maElVal = (maElVal + _elements.get(row).get(column));
    		_elements.get(row).put(column, maElVal);
    	}
    	else {
    		_elements.get(row).put(column, maElVal);
    	}
    	
    }
    
    public double getFactors() {
    	
    	//FOR TESTING ONLY, DO NOT USE
    	
    	return _operatorFactors;
    	
    }
    
    public Map<Integer, Double> getRow(int row) {
    	
    	//Returns the required row
    	
    	return _elements.get(row);
    	
    }
    
}
