/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mm;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import javax.imageio.ImageIO;
import javax.swing.JPanel;

/**
 *
 * @author manja
 */
public class Image extends JPanel{
    private BufferedImage shownImage = null;
    
    protected void paintComponent(Graphics g){
        super.paintComponent(g);
        if (shownImage != null){
            this.setPreferredSize(new Dimension(shownImage.getWidth(), shownImage.getHeight()));
            revalidate();
            g.drawImage(shownImage, 0, 0, null);
        }
    }
    
    public void openFile(File f) throws Exception{
        BufferedImage slika = ImageIO.read(f);
        shownImage = new BufferedImage(slika.getWidth(), slika.getHeight(), BufferedImage.TYPE_INT_RGB);
        for (int i=0; i<slika.getWidth(); i++){
            for (int j=0; j<slika.getHeight(); j++){
                shownImage.setRGB(i, j, slika.getRGB(i, j));
            }
        }
        repaint();
    }
}
