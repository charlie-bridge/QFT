package qft;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import utilities.Factorial;

public class InteractionMatrix {
    
    // Class for calculating and storing an interaction matrix associated with an interaction
    // term in the system Hamiltonian of general power.
    
    private int 		_systemSize;					
    private double 		_epsilon;						
    private double		_mass;							
    private int 		_numStates;						
    private int 		_phiPow;						// Power of the interaction term
    private int 		_rightMostP;					// Farthest right operator momentum
    private double 		_operatorFactors;				// Factors resulting from application of operators
    private int 		_momCombs;						// The number of momenta combinations
    private int 		_opCombs;						// The number of operator combinations
    private List<Map<Integer, Double>> 	_elements;		// Store of all the matrix elements
    
    public InteractionMatrix(int systemSize, double epsilon, double mass, int numStates, int phiPow) {
        
        // Constructor
        
        _systemSize = systemSize;
        _epsilon = epsilon;
        _mass = mass;
        _phiPow = phiPow;
        _numStates = numStates;
        _momCombs = (int)Math.pow(_systemSize, (phiPow-1));
        _opCombs = (int)Math.pow(2, phiPow);
        _elements = new ArrayList<Map<Integer, Double>>();
        
    }
    
    public void calcMatrix() {
        
        // Calculates the interaction matrix
        
        FockState rowState = new FockState(_systemSize, _epsilon, _mass);
        FockState opedState = new FockState(_systemSize, _epsilon, _mass);
        Integer column;
        int[] momenta = new int[_phiPow-1];        
        for(int row=0; row<_numStates; row++) {
            
            // Run through the matrix rows
        	// Set the row state and create the row Map
        	
            rowState.setAsIndex(row);
            Map<Integer, Double> rowMap = new HashMap<Integer, Double>();
    		_elements.add(row, rowMap);

        	// Reset the momenta combinations array
        	
        	for(int i=0; i<(_phiPow-1); i++) {
                momenta[i]=0;
            }            
            for(int j=0; j<_momCombs; j++) {
                
                // Run through the various momenta combinations
            
                for(int opType=0; opType<_opCombs; opType++) {
                    
                    // Run through the various operator combinations
                    
                    opedState.makeSameAs(rowState);
                    applyOperators(opedState, momenta, opType);                    
                    if(opedState.isValid()==true) {
                    	
                    	// Calculate the column
                    	
                    	column = opedState.calcIndex();                       
                    	if((column != null) && (column < _numStates)) {
                    		
                    		// If the resulting state is valid and within the cutoff, store the result
                    		
                    		store(row, column, momenta);
                        
                    	}                   	
                    }                
                }
            
                // Increment the momenta label array unless this is the final run through
                
                if(j != (_momCombs-1)) {
                	incrementMomentaLabel(momenta);
                }               
            }            
        }
        
    }
    
    public void incrementMomentaLabel(int[] momenta) {
    	
    	// PUBLIC FOR TESTING ONLY, DO NOT USE OUTSIDE THIS CLASS
        
        // Increments the momenta label array.
        // Add one to first value, but if that is = (_systemSize-1) then set it to zero
        // and increase the next one, but if that is = (_systemSize-1) etc.
        
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
        
        // If keepGoing is still true, all values must be N, so all combinations have been considered
        
    }
    
    public void applyOperators(FockState passedRow, int[] momenta, int opTypes) {
    	
    	// PUBLIC FOR TESTING ONLY, DO NOT USE OUTSIDE THIS CLASS
        
        // Applies operators to the passed row FockState according to the integers in momenta and opTypes.
        // We have 2^phiPow possibilities for the operators - so use binary:
        // for each 1 bit use annihilation, for each 0 bit use creation (on the passed <bra|)
        // (note: this means a 1 bit implies the operator is a creation operator).
        // The rightmost operator is special, its momentum is set by the others
        // this is stored outside this method for later use.
        
        _rightMostP = 0;
        _operatorFactors = 1.0;        
        for(int i=0; i<(_phiPow-1); i++) {
        	
        	// Run through the operators from left to right
            
            if(((opTypes & (1 << i)) != 0)) {
            	
            	// This bit is a 1, apply an annihilation operator to the relevant momentum mode
            	// and update the running operator factor and rightmost operator momentum.
            	// If this operation invalidates the state set the factor to zero.
                
                if(passedRow.getCoeff(momenta[i]) > 0) {
            		_operatorFactors = (_operatorFactors * Math.sqrt(passedRow.getCoeff(momenta[i])));
            		_operatorFactors = (_operatorFactors * Math.sqrt(2.0*FockState.calcFrequency(momenta[i], _systemSize, _epsilon, _mass)*_systemSize*_epsilon));                 
            	}
                else {
                	_operatorFactors = 0.0;
                }                
                passedRow.applyAnnihil(momenta[i]);
                
                // The rightmost momentum is set by the momentum just used
                
                if((momenta[i]%2) == 0) {
                	_rightMostP = _rightMostP - (int)Math.floor((momenta[i]+1)/2);               	
                }
                else {                	
                	_rightMostP = _rightMostP + (int)Math.floor((momenta[i]+1)/2);                	
                }                
            }
            else {
            	
            	// This bit is a 0, apply a creation operator to the relevant momentum mode
            	// and update the running operator factor and rightmost operator momentum.
                
            	_operatorFactors = (_operatorFactors * Math.sqrt(passedRow.getCoeff(momenta[i]) + 1.0));
            	_operatorFactors = (_operatorFactors * Math.sqrt(2.0*FockState.calcFrequency(momenta[i], _systemSize, _epsilon, _mass)*_systemSize*_epsilon));                      
                passedRow.applyCreation(momenta[i]);
                
                // The rightmost momentum is set by the momentum just used
                
                if((momenta[i]%2) == 0) {               	
                	_rightMostP = _rightMostP + (int)Math.floor((momenta[i]+1)/2);                	
                }
                else {                	
                	_rightMostP = _rightMostP - (int)Math.floor((momenta[i]+1)/2);                	
                }                
            }           
        }
        
        // Put the rightmost momentum in the required range before using it
        
        if(((opTypes & (1 << (_phiPow-1))) != 0)) {
        	
        	// This bit is a 1, apply an annihilation operator to the rightmost operator momentum mode.
        	// Complete the operator factor for this set of operators.
        	// If this operation invalidates the state set the factor to zero.
        	
        	while(_rightMostP >= (_systemSize/2)) {   			
    			_rightMostP = (_rightMostP - _systemSize);    				
    		}    		
    		while(_rightMostP < (-_systemSize/2)) {    				
    			_rightMostP = (_rightMostP + _systemSize);    				
    		}			
			if(_rightMostP >= 0) {				
				_rightMostP = (_rightMostP * 2);				
			}
			else {				
				_rightMostP = ((2*(-_rightMostP)) - 1);				
			}			
        	if(passedRow.getCoeff(_rightMostP) > 0) {        	    
        		_operatorFactors = (_operatorFactors * Math.sqrt(passedRow.getCoeff(_rightMostP)));
        		_operatorFactors = (_operatorFactors * Math.sqrt(2.0*FockState.calcFrequency(_rightMostP, _systemSize, _epsilon, _mass)*_systemSize*_epsilon));        		
        	}
        	else {
        		_operatorFactors = 0.0;
        	}        	
        	passedRow.applyAnnihil(_rightMostP);            
        }
        else {
        	
        	// This bit is a 0, apply a creation operator to the rightmost operator momentum mode.
        	// Complete the operator factor for this set of operators.
        	
        	_rightMostP = (_rightMostP * -1);       	
        	while(_rightMostP >= (_systemSize/2)) {   			
    			_rightMostP = (_rightMostP - _systemSize);    				
    		}    		
    		while(_rightMostP < (-_systemSize/2)) {   				
    			_rightMostP = (_rightMostP + _systemSize);    				
    		}       	
        	if(_rightMostP >= 0) {				
				_rightMostP = (_rightMostP * 2);				
			}
			else {				
				_rightMostP = ((2*(-_rightMostP)) - 1);				
			}        	
        	_operatorFactors = (_operatorFactors * Math.sqrt(passedRow.getCoeff(_rightMostP) + 1.0));
        	_operatorFactors = (_operatorFactors * Math.sqrt(2.0*FockState.calcFrequency(_rightMostP, _systemSize, _epsilon, _mass)*_systemSize*_epsilon));       	
            passedRow.applyCreation(_rightMostP);               
        }
        
    }
    
    private void store(int row, int column, int[] momenta) {
        
        // Stores the passed matrix value.
    	// First, calculate the value to store.
    	// Note: The rightmost momentum value has already been calculated by the applyOperators method.
        // Because calcMatrix goes through the rows in order, they will be created in order.
    	
    	double maElVal = (Math.pow((_systemSize*_epsilon), -(_phiPow-1)) * Math.pow(2.0, -_phiPow) * (1.0/Factorial.calc(_phiPow)));
    	maElVal = (maElVal * _operatorFactors);
    	maElVal = (maElVal / FockState.calcFrequency(_rightMostP, _systemSize, _epsilon, _mass));    	
    	for(int i=0; i<(_phiPow-1); i++) {
    		
    		maElVal = (maElVal / FockState.calcFrequency(momenta[i], _systemSize, _epsilon, _mass));   		
    	}
    	
    	// Check if a value already exists at this location, if it does, just add to it
    	
    	if(_elements.get(row).containsKey(column) == true) {
    		maElVal = (maElVal + _elements.get(row).get(column));
    		_elements.get(row).put(column, maElVal);	
    	}
    	else {
    		_elements.get(row).put(column, maElVal);
    	}
    	
    }
    
    public double getFactors() {
    	
    	// FOR TESTING ONLY, DO NOT USE
    	
    	return _operatorFactors;
    	
    }
    
    public Map<Integer, Double> getRow(int row) {
    	
    	// Returns the required row
    	
    	return _elements.get(row);
    	
    }
    
}
