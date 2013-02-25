package qft;

public class IntMatPos {
    
    //Class to use as key for Matrix map
    
    private final int _row;
    private final int _column;
    
    public IntMatPos(int rowValue, int columnValue) {
        
        _row = rowValue;
        _column = columnValue;
        
    }
    
    public int getRow() {
        
        return _row;
        
    }
    
    public int getColumn() {
        
        return _column;
        
    }
    
    public int hashCode() {
        return _row ^ _column;
    }
    
    public boolean equals(Object objectPassed) { 
      IntMatPos keyPassed = (IntMatPos)objectPassed;
      return objectPassed instanceof IntMatPos && keyPassed.getRow() == _row && keyPassed.getColumn() == _column;

    }
    
}
