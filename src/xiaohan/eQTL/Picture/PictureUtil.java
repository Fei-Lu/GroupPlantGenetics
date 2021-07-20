package xiaohan.eQTL.Picture;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.image.BufferedImage;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import com.sun.image.codec.jpeg.JPEGCodec;
import com.sun.image.codec.jpeg.JPEGImageEncoder;

public class PictureUtil extends JFrame implements MouseListener {

    public PictureUtil() {
        setSize(1600, 1200);// 根据要求调整大小
        setLocation(200, 100);
        setTitle("获得图片上任意点坐标");
        setResizable(false);

        Container con = getContentPane();

        ImageIcon bgIcon = new ImageIcon("/Users/yxh/Documents/eQTL/SiPAS/RIN/7.6.jpeg");// 注意图片的路径
        ImagePanel backpicPanel = new ImagePanel(bgIcon,2);
        backpicPanel.addMouseListener(this);
        con.add(backpicPanel, BorderLayout.CENTER);

        tipLabel = new JLabel(
                "--------------------提示：坐标依次打印在屏幕上!--------------------");
        con.add(tipLabel, BorderLayout.SOUTH);
    }

    private JLabel tipLabel;

    private static void createImage(String fileLocation, BufferedImage image) {
        try {
            FileOutputStream fos = new FileOutputStream(fileLocation);
            BufferedOutputStream bos = new BufferedOutputStream(fos);
            JPEGImageEncoder encoder = JPEGCodec.createJPEGEncoder(bos);
            encoder.encode(image);
            bos.close();
            fos.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        PictureUtil frame = new PictureUtil();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }

    public void mouseClicked(MouseEvent arg0) {
        // TODO Auto-generated method stub

    }

    public void mouseEntered(MouseEvent arg0) {
        // TODO Auto-generated method stub

    }

    public void mouseExited(MouseEvent arg0) {
        // TODO Auto-generated method stub

    }

    public void mousePressed(MouseEvent e) {
        // TODO Auto-generated method stub
        int x = e.getX();
        int y = e.getY();
        String message = "(" + x + "," + y + ")";
        tipLabel.setText(message);
        System.out.println(message);
    }

    public void mouseReleased(MouseEvent arg0) {
        // TODO Auto-generated method stub

    }
}

class ImagePanel extends JPanel {
    private Image img;

    public ImagePanel(ImageIcon imageIcon,double ratio) {
        int width = (int)(imageIcon.getIconWidth()*ratio);
        int height = (int)(imageIcon.getIconHeight()*ratio);
        img = imageIcon.getImage().getScaledInstance(width, height, Image.SCALE_DEFAULT);
    }

    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        g.drawImage(img, 0, 0, this);
    }

}