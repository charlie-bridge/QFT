package utilities;

public class Factorial {
    
    // Calculates the factorial of a number
    // Not particularly hardy, intended for low numbers, use with care.
	
	private int _value;
	
	public Factorial(int n) {
		
		_value = calc(n);
		
	}
	
	public int getValue() {
		
		// Gets the value
		
		return _value;
		
	}
	
	public static int calc(int n) {
		
		// Calculates the value of n factorial
		
		int factorial = 1;
		for(int i=1; i<=n; i++) {			
			factorial = (factorial * i);			
		}
		
		return factorial;
		
	}

}
