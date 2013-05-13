package qft;

public class FockState {
    
    // Class for one of the members of the infinite Fock Space
    // Uses Cantor's Pairing Function to map a state to an index (and vice versa)	
	// Even coefficients refer to forward moving particles,
	// odd coefficients refer to backward moving particles
    
    private Integer[] _coeffs;
    private int _systemSize;
    private double _epsilon;
    private double _mass;
    private boolean _isValid;
    
    public FockState(int systemSize, double epsilon, double mass) {
        
        // Constructor, initialise to the vacuum state
        
        _systemSize = systemSize;
        _epsilon = epsilon;
        _mass = mass;
        _coeffs = new Integer[systemSize];
        setAsVacuum();
        _isValid = true;
        
    }
    
    public void setAsVacuum() {
        
        // Sets all coefficients to zero - the Free Hamiltionian vacuum Fock State
        
        for(int i=0; i<_systemSize; i++) {
            _coeffs[i] = 0;
        }
        
    }
    
    public Integer getCoeff(int x) {
        
        // Gets the xth component of the state, returns null if x is outside the
        // required range
        
        if((x<_systemSize) && (x>=0)) {
            return _coeffs[x];
        }
        else {
            return null;
        }
        
    }

    public void setCoeff(int y, int z) {
        
        // Sets y component of the state to z, does nothing if y is outside the 
        // required range or z is negative
        
        if( (y<_systemSize) && (y>=0) && (z>=0)) {
            _coeffs[y] = z;
        }
        
    }
    
    public Integer calcIndex() {
        
        // Calculates and returns the index for this Fock State
        
        long index = 0;
        long a;
        long b;        
        for(int i=0; i<(_systemSize-1); i++) {
            
            // The index is calculated by using Cantor's Pairing Functions recursively,
            // see p32-34 of lab book
            
            a = _coeffs[_systemSize - 2 - i];
            if(i==0) {
                b = _coeffs[_systemSize - 1];
            }
            else {
                b = index;
            }            
            index = ((((a+b)*(a+b+1))/2) + b);            
            if(index > Integer.MAX_VALUE)                
            	return null;      
        }
        
        return (int)index;
        
    }
    
    public void setAsIndex(int index) {
        
        // Sets this Fock State to the one which matches the passed index
        
        setAsVacuum();
        int rollingIndex = index;
        long w;        
        for(int i=0; i<(_systemSize-1); i++) {
            
            // The Fock State is calculated using Cantor's Pairing Functions,
            // see p32-34 of lab book
            
            w = (long)Math.floor((Math.sqrt((8*(double)rollingIndex)+1.0) - 1)/2);           
            rollingIndex = (int)(rollingIndex - ((w*(w+1))/2));
            _coeffs[i] = (int)(w - rollingIndex);            
        }        
        _coeffs[_systemSize - 1] = (int)rollingIndex;
        
    }
    
    public void applyCreation(int p) {
        
        // Applies a creation operator for the pth momentum mode
        // For even p this creates a forward moving particle
    	// For odd p this creates a backwards moving particle
        
        if((p<_systemSize) && (p>=0)) {
            _coeffs[p]++;
        }
        
    }
    
    public void applyAnnihil(int r) {
        
        // Applies an annihilation operator for the rth momentum mode,
        // does nothing if r is outside the required range.
        // If the annihilation operator destroys the state, set it to invalid.   	
    	// For even p this destroys a forward moving particle
    	// For odd p this destroys a backwards moving particle
        
        if((r<_systemSize) && (r>=0)) {        	
            _coeffs[r]--;          
            if(_coeffs[r]==-1) {
                _isValid = false;
            }           
        }
                
    }
    
    public boolean isValid() {
        
        // Returns true if the state is valid, otherwise returns false
        
        return _isValid;
        
    }
    
    public int getSystemSize() {
        
        // Gets the system size
        
        return _systemSize;
        
    }
    
    public void makeSameAs(FockState passedState) {
        
        // Makes this Fock state the same as the passed one    	
    	// Note: does not check mass or epsilon
        
        if(passedState.getSystemSize() == _systemSize) {
            for(int i=0; i<_systemSize; i++) {
                _coeffs[i] = passedState.getCoeff(i);
            }            
            _isValid = passedState.isValid();            
        }
        
    }
    
    public double calcEnergy() {
    	
    	// Returns the energy of this Fock State
    	
    	double energy = 0.0;    	
    	for(int i=0; i<_systemSize; i++) {    		
    		if(_coeffs[i] > 0) {
    			energy = (energy + (_coeffs[i]*calcFrequency(i, _systemSize, _epsilon, _mass)));
    		}    		
    	}
    	
    	return energy;
    	
    }
    
    public static double calcFrequency(int p, int systemSize, double epsilon, double mass) {
    	
    	// Returns the frequency of a particular momentum mode.    	
    	// First, find the appropriate value to put in the formula (see p43 of lab book).
    	
    	int mode = (int)Math.floor((p+1)/2);
    	
    	return (Math.sqrt(((4/(Math.pow(epsilon, 2))) * Math.pow(Math.sin((Math.PI*(double)mode)/((double)systemSize)),  2)) + Math.pow(mass, 2)));
    	
    }
    
    public static Integer getIndex1PState(int p, int systemSize, double epsilon, double mass) {
    	
        // Gets the index of the one-particle state with a particle at 
        // the index p as required for the Fourier transform.
    	// Epsilon and mass are not used in this class other than to construct the state.   	
    	// Limit of p for forward moving particles calculated by rounding down (sysSize-1)/2.
    	// Minus one comes from starting index at zero, rounding down due to backwards particles
    	// being indexed before their forward counterparts.
    	
    	FockState state = new FockState(systemSize, epsilon, mass);
    	
    	//Check p is in the right range, if not return null
    	
    	if((p>=0) && (p<systemSize)) {    	
    		if(p <= (int)Math.floor((systemSize-1)/2)) {
    		
    			// For p up to the limit create a forward moving particle
    		
    			state.setCoeff((2*p), 1);   		
    		}
    		else {
    		
    			// For p beyond the limit create a backward moving particle
    			
    			state.setCoeff(((2*(systemSize-p))-1), 1);   		
    		}    	
    	}
    	else {    		
    		return null;   		
    	}
    	
    	return state.calcIndex();
    	
    }
    
    public static Integer getIndex2PState(int p, int q, int systemSize, double epsilon, double mass) {
    	
    	// Gets the index of the two-particle state with particles at the 
    	// indexes p and q as required for the Fourier transform.
        // Epsilon and mass are not used in this class other than to construct the state.
    	
    	FockState state = new FockState(systemSize, epsilon, mass);
    	
    	// Perform the same analysis as the one-particle state twice.   	
    	// Check p and q are in the right range, if not return null.
    	
    	if((p>=0) && (p<systemSize) && (q>=0) && (q<systemSize)) {    	
    		if(p <= (int)Math.floor((systemSize-1)/2)) {
    		
    			// For p up to the limit create a forward moving particle
    		
    			state.applyCreation(2*p);    		
    		}
    		else {
    		
    			// For p beyond the limit create a backward moving particle
    			
    			state.applyCreation((2*(systemSize-p))-1);   		
    		}    	    	
        	if(q <= (int)Math.floor((systemSize-1)/2)) {
        		
        		// For q up to the limit create a forward moving particle
        		
        		state.applyCreation(2*q);        		
        	}
        	else {
        		
        		// For q beyond the limit create a backward moving particle

        		state.applyCreation((2*(systemSize-q))-1);        		
        	}    	
    	}
    	else {    		
    		return null;    		
    	}
    	
    	return state.calcIndex();
    	
    }
    
}
