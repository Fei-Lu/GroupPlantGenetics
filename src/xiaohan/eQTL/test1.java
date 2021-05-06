import javax.swing.*;
import java.awt.*;
public class test1 {
    //X- and Y-Coordinates
    public int x, y;

    //Player graphic
    Image img = Toolkit.getDefaultToolkit().getImage("path.to.your.image.file");

    public test1() {

    }

    //Rendering method
    //You should call this in your frame class where you override paint(Graphics g)
    //or paintComponent(Graphics g)
    //
    //like this: player.draw(g,this)
    public void draw(Graphics g, JPanel jp) {
        //Draw the Image at the right position
        g.drawImage(img, x, y, jp);
    }

    //Get X coordinate
    public int getX() {
        return x;
    }

    //Get Y coordinate
    public int getY() {
        return y;
    }

    //Update method
    //You should call this in your game loop
    public void update() {
        //Handle button presses, collision testing or other stuff relating
        //the player object here
    }
}