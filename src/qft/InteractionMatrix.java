package qft;

import qft.IntMatPos;

import java.util.*;

public class InteractionMatrix {
    
    //Class for calculating and storing an interaction matrix for theory of general interaction power
    
    private int _systemSize;
    private int _cutOffMom;
    private int _numStates;
    private int _phiPow;
    private HashMap<IntMatPos, Double> _elements;
    
    public InteractionMatrix(int systemSize, int cutOffMom, int phiPow) {
        
        //Cut off momentum should be greater than or equal to one
        
        _systemSize = systemSize;
        _cutOffMom = cutOffMom;
        _phiPow = phiPow;
        _numStates = (((cutOffMom*(cutOffMom+1))/2) + 1);
        _elements = new HashMap<IntMatPos, Double>(_numStates*(int)Math.pow(systemSize, 2)*(int)Math.pow(2, (phiPow-1)));
        
    }
    
    public void calcMatrix() {
        
        FockState rowState = new FockState(_systemSize);
        FockState opedState = new FockState(_systemSize);
        int column;
        
        int[] momenta = new int[_phiPow-1];
        for(int i=0; i<_phiPow-1; i++) {
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
        
        int rightMostP = 0;
        
        for(int i=0; i<(_phiPow-2); i++) {
            
            if(((opTypes & (1 << i)) != 0)) {
                
                passedRow.applyCreation(momenta[i]);
                rightMostP = rightMostP + momenta[i];
                
            }
            else {
                
                passedRow.applyAnnihil(momenta[i]);
                rightMostP = rightMostP - momenta[i];
                
            }
            
        }
        
        //Put the rightmost momenta in the required range before using it
        
        if(((opTypes & (1 << (_phiPow-1))) != 0)) {
            
            rightMostP = Math.abs(rightMostP) % _systemSize;
            passedRow.applyCreation(rightMostP);
            
        }
        else {
            
            rightMostP = _systemSize - (Math.abs(rightMostP) % _systemSize);
            passedRow.applyAnnihil(rightMostP);
            
        }
        
    }
    
    private void store(int row, int column, int[] momenta) {
        
        //Stores the matrix element, check for already existing elements with the same position
        
        
    }
    

}
