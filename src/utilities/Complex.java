package utilities;

public class Complex {
    
    // SUPERCEDED BY CODEBASE VERSION OF COMPLEX CLASS, DO NOT USE
	
	private double _real;
	private double _comp;
	
	public Complex(double real, double comp) {
		_real = real;
		_comp = comp;
	}
	
	public double realPart() {
		return _real;
	}
	
	public double complexPart() {
		return _comp;
	}
	
	public void setReal(double inputReal) {
		_real = inputReal;
	}
	
	public void setComp(double inputComp) {
		_comp = inputComp;
	}
	
	public void add(Complex addor) {
		_real += addor._real;
		_comp += addor._comp;
	}
	
	public void subtract(Complex subtractor) {
		_real -= subtractor._real;
		_comp -= subtractor._comp;
	}
	
	public double modulus() {
		return Math.sqrt(Math.pow(_real, 2) + Math.pow(_comp, 2));
	}
	
}
