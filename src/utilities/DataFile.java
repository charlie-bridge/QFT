package utilities;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class DataFile {
    
    // Gets a file ready for writing to. Sets the headings as desired.
    // File output intended for testing.
    
    private File _file;
    
    public DataFile(String fileSuffix) throws IOException {
        
        // Constructor with file name
        
        _file = new File("C:\\Users\\Charlie\\Desktop\\QFTFiles\\" + fileSuffix + ".log");        
        if (!_file.exists()) {
            _file.createNewFile();
        }
        
    }
    
    public void writeHeadings(String[] headings) throws IOException {
        
      // Writes a line of headings, don't make them really long
   
      BufferedWriter line = new BufferedWriter(new FileWriter(_file.getAbsoluteFile(), true));            
      for(int i=0; i<headings.length; i++) {                
          line.append(headings[i]);               
              if(i < (headings.length - 1)) {
                  line.append("\t");
              }                
      }            
      line.close();
        
    }
    
    public void writeLine(Double[] data) throws IOException {
        
        // Writes a line of data, watch out for long doubles

        BufferedWriter line = new BufferedWriter(new FileWriter(_file.getAbsoluteFile(), true));      
        line.newLine();
        String dataString;            
        for(int i=0; i<data.length; i++) {                
            dataString = Double.toString(data[i]);
            line.append(dataString);                
            if(i < (data.length - 1)) {
                line.append("\t");
            }                
        }            
        line.close();
        
    }

}
