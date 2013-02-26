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
        for(int i=0; i<(_phiPow-1); i++) {
            momenta[i]=0;
        }
        
        for(int row=0; row<_numStates; row++) {
            
            //Run through the matrix rows
            
            rowState.setAsIndex(row);
            
            for(int j=0; j<((int)Math.pow(2, (_phiPow-1))); j++) {
                
                //Run through the various momenta combinations
            
                for(int k=0; k<(_phiPow-1); k++) {
                    
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
    
    private void incrementMomentaLabel(int[] momenta) {
        
        //Increments the momenta label array
        //Add one to first value, but if that is = (_systemSize-1) then set it to zero
        //increase the next one, but if that is = (_systemSize-1) etc.
        
        boolean keepGoing = true;
        
        for(int i=0; (i<(_phiPow-1) && (keepGoing == true)); i++) {
            
            if(momenta[i]<_systemSize) {
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
    
    private void applyOperators(FockState passedRow, int[] momenta, int opTypes) {
        
        //Applies operators to the passed row FockState according to the integer i
        //We have 2^(phiPow-1) possibilities for the operators - so use binary
        //For each 0 bit use annihilation, for each 1 bit use creation
        //The rightmost operator is special, its momentum is set by the others
    	//this is stored outside this method for use later
        
        _rightMostP = 0;
        _operatorFactors = 1.0;
        
        for(int i=0; i<(_phiPow-2); i++) {
            
            if(((opTypes & (1 << i)) != 0)) {
                
            	_operatorFactors = (_operatorFactors * Math.sqrt(passedRow.getCoeff(momenta[i]) + 1.0));
                passedRow.applyCreation(momenta[i]);
                _rightMostP = _rightMostP + momenta[i];
                
            }
            else {
                
            	if(passedRow.getCoeff(momenta[i]) > 0) {
            		_operatorFactors = (_operatorFactors * Math.sqrt(passedRow.getCoeff(momenta[i])));
            	}
                passedRow.applyAnnihil(momenta[i]);
                _rightMostP = _rightMostP - momenta[i];
                
            }
            
        }
        
        //Put the rightmost momenta in the required range before using it
        
        if(((opTypes & (1 << (_phiPow-1))) != 0)) {
            
            _rightMostP = Math.abs(_rightMostP % _systemSize);
            passedRow.applyCreation(_rightMostP);
            
        }
        else {
            
            _rightMostP = _systemSize - (Math.abs(_rightMostP % _systemSize));
            passedRow.applyAnnihil(_rightMostP);
            
        }
        
    }
    
    private void store(int row, int column, int[] momenta) {
        
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
    

}
