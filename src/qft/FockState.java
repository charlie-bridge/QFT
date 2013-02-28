package qft;

public class FockState {
    
    //Class for one of the members of the infinite Fock Space
    //includes two methods of indexing: one for quickly finding states or creating them,
    //one for ordering them in a way that is useful for Fourier analysis
    
    private int[] _coeffs;
    private int _systemSize;
    private double _epsilon;
    private double _mass;
    private boolean _isValid;
    
    public FockState(int systemSize, double epsilon, double mass) {
        
        //Initialise the Fock State, set to the vacuum state
        
        _systemSize = systemSize;
        _epsilon = epsilon;
        _mass = mass;
        _coeffs = new int[systemSize];
        setAsVacuum();
        _isValid = true;
        
    }
    
    public void setAsVacuum() {
        
        //Sets all coefficients to zero - the vacuum Fock State
        
        for(int i=0; i<_systemSize; i++) {
            _coeffs[i] = 0;
        }
        
    }
    
    public int getCoeff(int x) {
        
        //Gets the xth component of the state,
        //returns -1 if x is outside the required range
        
        if((x<_systemSize) && (x>=0)) {
            return _coeffs[x];
        }
        else {
            return -1;
        }
        
    }

    public void setCoeff(int y, int z) {
        
        //Sets y component of the state to z,
        //does nothing if y is outside the required range or z is negative
        
        if( (y<_systemSize) && (y>=0) && (z>=0)) {
            _coeffs[y] = z;
        }
        
    }
    
    public int calcIndex() {
        
        //Calculates and returns the index for this Fock State
        
        int index = 0;
        int a;
        int b;
        
        for(int i=0; i<(_systemSize-1); i++) {
            
            //The index is calculated by using Cantor's Pairing Functions recursively, see p32-34
            
            a = _coeffs[_systemSize - 2 - i];
            if(i==0) {
                b = _coeffs[_systemSize - 1];
            }
            else {
                b = index;
            }
            
            index = ((((a+b)*(a+b+1))/2) + b);
            
        }
        
        return index;
        
    }
    
    public void setAsIndex(int index) {
        
        //Sets this Fock State to the one which matches the passed index
        
        setAsVacuum();
        int rollingIndex = index;
        long w;
        
        for(int i=0; i<(_systemSize-1); i++) {
            
            //The Fock State is calculated using Cantor's Pairing Functions, see p32-34
            
            w = (long)Math.floor((Math.sqrt((8*(double)rollingIndex)+1.0) - 1)/2);
            
            rollingIndex = (int)(rollingIndex - ((w*(w+1))/2));
            _coeffs[i] = (int)(w - rollingIndex);
            
        }
        
        _coeffs[_systemSize - 1] = (int)rollingIndex;
        
    }
    
    public void applyCreation(int p) {
        
        //Applies a creation operator for the pth momentum mode
        //does nothing if p is outside the required range
        
        if((p<_systemSize) && (p>=0)) {
            _coeffs[p]++;
        }
        
    }
    
    public void applyAnnihil(int r) {
        
        //Applies an annihilation operator for the rth momentum mode
        //does nothing if r is outside the required range
        //If the annihilation operator destroys the state, set it to invalid
        
        if((r<_systemSize) && (r>=0)) {
        	
            _coeffs[r]--;
            
            if(_coeffs[r]==-1) {
                _isValid = false;
            }
            
        }
                
    }
    
    public boolean isValid() {
        
        //Returns true if the state is valid, otherwise returns false
        
        return _isValid;
        
    }
        
    public void incrementState() {
        
        //Increments the Fock State according to the rules laid out on p17
        //This method is used to supply Fourier transforms with the coefficients in the correct order
        
        //First, find the last and second last non-zero numbers in the state
        
        int lastNumPos = 0;
        int lastNumVal = _coeffs[0];
        int scndLastNumPos = 0;
        int scndLastNumVal = 0;
        
        for(int k=1; k<_systemSize; k++){
            
            if(_coeffs[k] != 0) {
                scndLastNumPos = lastNumPos;
                scndLastNumVal = lastNumVal;
                lastNumPos = k;
                lastNumVal = _coeffs[k];
            }
            
        }
        
        if(lastNumPos != (_systemSize - 1)) {
            
            //If the last number is not at the end of the state, take one
            //from it and place that one just ahead of it
            
            _coeffs[lastNumPos] = (lastNumVal - 1);
            _coeffs[lastNumPos + 1] = 1;
            
        }
        else {
            
            if(scndLastNumVal != 0) {
                
                //If the last number is at the end of the state (and we are not
                //in the special case where that is the only number) then take
                //one off the 2nd last number, add it to the last and place it
                //just ahead of the 2nd last
                
                _coeffs[scndLastNumPos] = (scndLastNumVal - 1);
                _coeffs[lastNumPos] = 0;
                _coeffs[scndLastNumPos + 1] = (lastNumVal + 1);
                
            }
            else {
                
                //In this special case, move to the start of the next set of
                //state which includes one more particle
                _coeffs[lastNumPos] = 0;
                _coeffs[0] = lastNumVal + 1;
                
            }
            
        }
        
    }

    public int getSystemSize() {
        
        //Gets the system size
        
        return _systemSize;
        
    }
    
    public void makeSameAs(FockState passedState) {
        
        //Makes this Fock state the same as the passed one
    	
    	//DOES NOT CHECK MASS OR EPSILON
        
        if(passedState.getSystemSize() == _systemSize) {
            for(int i=0; i<_systemSize; i++) {
                _coeffs[i] = passedState.getCoeff(i);
            }
            
            _isValid = passedState.isValid();
            
        }
        
    }
    
    public double calcEnergy() {
    	
    	//Returns the energy of this Fock State
    	
    	double energy = 0;
    	
    	for(int i=0; i<_systemSize; i++) {
    		
    		if(_coeffs[i] > 0) {
    			energy = (energy + (_coeffs[i]*calcFrequency(i, _systemSize, _epsilon, _mass)));
    		}
    		
    	}
    	
    	return energy;
    	
    }
    
    public static double calcFrequency(int p, int systemSize, double epsilon, double mass) {
    	
    	//Returns the frequency of a particular momentum mode
    	
    	return (Math.sqrt(((4/(Math.pow(epsilon, 2))) * Math.pow(Math.sin((Math.PI*(double)p)/((double)systemSize)),  2)) + Math.pow(mass, 2)));
    	
    }
    
}
