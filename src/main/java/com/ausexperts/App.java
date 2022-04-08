package com.ausexperts;
import com.ausexperts.gas.*;
/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args )
    {
        System.out.println( "Hello World!" );
        Gas g = new Gas(0.72,0,0);
        double x = g.getGaszfactor(2000, 140);
        System.out.println();
    }
}
