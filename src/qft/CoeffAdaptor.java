package qft;

import uk.ac.cam.cal56.maths.Complex;

public class CoeffAdaptor {

    public CoeffAdaptor() {
        
        //Constructor does nothing
        
    }
    
    public static Complex[] setCoeffs(int systemSize, int numStates, Complex[] naiveCoeffs) {
        
        //Returns the reordered coefficients
        
        Complex[] newCoeffs = new Complex[numStates];
        for(int i=0; i<numStates; i++) {
            newCoeffs[i] = Complex.zero();
        }
        
        //Set the vacuum state coefficient
        
        newCoeffs[0] = naiveCoeffs[0];
        
        //Set the one-particle coefficients
        
        FockState fockState = new FockState(systemSize, 1, 1);
        Integer index;
        
        for(int j=0; j<systemSize; j++) {
            
            if(j<(systemSize/2)) {
                fockState.applyCreation(2*j);
            }
            else {
                fockState.applyCreation((2*(systemSize-j))-1);
            }
            
            index = fockState.calcIndex();
            
            if((index != null) && (index<numStates)) {
                newCoeffs[index] = new Complex(naiveCoeffs[j+1].real(), naiveCoeffs[j+1].imag());
            }
            
            fockState.setAsVacuum();
            
        }
        
        //Set the two-particle coefficients
        //Can we do this within the single particle cycle??
        
        int naiveIndex = (systemSize + 1);
        
        for(int k=0; k<systemSize; k++) {
            
            for(int h=k; h<systemSize; h++){

                if(k<(systemSize/2)) {
                    fockState.applyCreation(2*k);
                }
                else {
                    fockState.applyCreation((2*(systemSize-k))-1);
                }
                
                if(h<(systemSize/2)) {
                    fockState.applyCreation(2*h);
                }
                else {
                    fockState.applyCreation((2*(systemSize-h))-1);
                }
                
                index = fockState.calcIndex();
                
                if((index != null) && (index<numStates)) {
                    newCoeffs[index] = new Complex(naiveCoeffs[naiveIndex].real(), naiveCoeffs[naiveIndex].imag());
                }
                
                naiveIndex++;
                
                fockState.setAsVacuum();
                
            }
            
        } 
        
        return newCoeffs;
        
    }
    
}
