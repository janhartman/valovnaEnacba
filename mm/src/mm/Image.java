/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mm;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import javax.imageio.ImageIO;
import java.awt.color.ColorSpace;
import javax.swing.JPanel;

/**
 *
 * @author manja
 */
public class Image extends JPanel{
    private BufferedImage shownImage = null;
    
    private final static double MIN = 0 ;
    private final static double MAX = 5 ;

    
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
    
    private int getColorForValue(float value) {
        float BLUE_HUE = Color.RGBtoHSB(0, 0, 255, null)[0];
        float RED_HUE = Color.RGBtoHSB(255, 0, 0, null)[0];
        
        if (value < MIN || value > MAX) {
            return 0 ;
        }
        float hue = (float) (BLUE_HUE + (RED_HUE - BLUE_HUE) * (value - MIN) / (MAX - MIN) );
        return Color.HSBtoRGB(hue, (float) 1.0, (float) 1.0);
    }
    
    
    public void colorImage(int h, int w, float[][] matrix) throws Exception{
        int[][] colors = new int[matrix.length][matrix.length];
        //matrix = null;
        
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                colors[i][j] = getColorForValue(matrix[i][j]);
            }
        }
        
        shownImage = null;
        if(matrix == null){
            shownImage = new BufferedImage(matrix.length, matrix.length, BufferedImage.TYPE_INT_RGB);
            for (int i=0; i<w; i++){
                for (int j=0; j<h; j++){
                    shownImage.setRGB(i, j,(i*j)%255);
                }
            }
        }
        else{
            shownImage = new BufferedImage(w,h, BufferedImage.TYPE_INT_RGB);
            for (int i=0; i< matrix.length; i++){
                for (int j=0; j< matrix.length; j++){
                    shownImage.setRGB(i, j, colors[i][j]);
                }
            }
        }
        
        
        repaint();
    }
}
