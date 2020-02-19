/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package feilu;

/**
 *
 * @author feilu
 */
public class Start {
    
    public Start () {
//        this.testCallableFuture();
//        this.bitSet();
        this.listAndArray();
    }

    public void listAndArray() {
        new ListAndArray();
    }

    public void bitSet () {
        new BitSetTest();
    }

    public void testCallableFuture () {
        new CallableFuture();
    }
     
    public static void main (String[] args) {
        new Start();
    }
    
}


