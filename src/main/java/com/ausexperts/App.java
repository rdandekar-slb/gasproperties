package com.ausexperts;

import com.ausexperts.gas.Gas;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args )
    {
        //System.out.println( "Hello World!" );
        Gas g = new Gas(0.7);
        System.out.println(g.getSpecificGravity());
    }
}
