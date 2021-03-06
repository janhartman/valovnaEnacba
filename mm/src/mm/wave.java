/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mm;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author manja
 */
public class wave extends javax.swing.JFrame {

    /**
     * Creates new form wave
     */
    public wave() {
        initComponents();
        //velikost mreze
        n = 30;

        //pozicije
        u = new double[n][n];

        //hitrosti
        v = new double[n][n];

        //utezi
        w = new double[n][n];

        //zacetni pogoji za u in v
        //definiraj custom
        v[n / 2][n / 2] = - 0.1;
        u[n/2][n/2] = 1;
        //v[2*n/3][2*n/3] = 0.5;
        

        //hitrost valovanja
        c =  0.1;

        //faktor dusenja
        k =  0.001;

        //razlike med tockami
        h =  0.2;

        //pogoj za rob
        edgeCondition =  0.0;
        
        //ali animacija tece
        animating = false;
        
        currentTime = System.currentTimeMillis();
    }

    //velikost mreze
    int n;

    //pozicije, hitrosti, utezi
    double[][] u, v, w;
    
    //parametri
    double c, k, h, edgeCondition;
    
    boolean animating;
    
    long currentTime;
    
    
    //metoda za racunanje novih hitrosti in pozicij
    public void update() {

        double[][] newU = new double[n][n];
        double[][] newV = new double[n][n];
        double[][] a = new double[n][n];
        long newTime = System.currentTimeMillis();
        int dt = (int) (newTime - currentTime);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {

                //pogoji za rob
                //nastavljamo
                if (i == 0 || j == 0 || i == n - 1 || j == n - 1) {
                    a[i][j] = edgeCondition;
                    newU[i][j] = edgeCondition;
                    newV[i][j] = edgeCondition;
                } 

                //sicer
                else {
                    a[i][j] = c * c / (h * h) * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - 4 * u[i][j]) - ( k  * v[i][j]);

                    
                    //nova hitrost je stara + sprememba
                    newV[i][j] = v[i][j] + a[i][j] * (double)dt / 1000.0 ;

                    //nova pozicija je stara + sprememba
                    newU[i][j] = u[i][j] + v[i][j] * (double)dt / 1000.0;
                }
            }
        }
        
        currentTime = newTime;
        
        

        u = newU;
        v = newV;
        
      // printMatrix();
        
        
        
        try{
            image.colorImage(image.getHeight(), image.getWidth(), u);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void runAnimation() {
        new Thread(new Runnable() {
            public void run() {
                
                int i = 100;
                while(true) {
                    update();
                    
                    try{
                        Thread.sleep(50);
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                    }
                 }
              }
        }).start();
        
    }
    
    public void printMatrix () {
         
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.printf("%.2f ", u[i][j]);
            }
            System.out.println();
        }
        
        System.out.println();
        
    }
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jPanel1 = new javax.swing.JPanel();
        imagePanel = new javax.swing.JPanel();
        image = new mm.Image();
        jLabel6 = new javax.swing.JLabel();
        controls = new javax.swing.JPanel();
        jLabel2 = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        input_constant = new javax.swing.JTextField();
        input_k = new javax.swing.JTextField();
        input_speed = new javax.swing.JTextField();
        jLabel5 = new javax.swing.JLabel();
        buttonOK = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("Wave equation visualization");

        image.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(102, 102, 102)));
        image.addComponentListener(new java.awt.event.ComponentAdapter() {
            public void componentResized(java.awt.event.ComponentEvent evt) {
                imageComponentResized(evt);
            }
        });

        javax.swing.GroupLayout imageLayout = new javax.swing.GroupLayout(image);
        image.setLayout(imageLayout);
        imageLayout.setHorizontalGroup(
            imageLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        imageLayout.setVerticalGroup(
            imageLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        jLabel6.setText("Status bar");

        javax.swing.GroupLayout imagePanelLayout = new javax.swing.GroupLayout(imagePanel);
        imagePanel.setLayout(imagePanelLayout);
        imagePanelLayout.setHorizontalGroup(
            imagePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(imagePanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel6, javax.swing.GroupLayout.DEFAULT_SIZE, 468, Short.MAX_VALUE))
            .addGroup(imagePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, imagePanelLayout.createSequentialGroup()
                    .addContainerGap()
                    .addComponent(image, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
        );
        imagePanelLayout.setVerticalGroup(
            imagePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, imagePanelLayout.createSequentialGroup()
                .addGap(0, 0, Short.MAX_VALUE)
                .addComponent(jLabel6, javax.swing.GroupLayout.PREFERRED_SIZE, 26, javax.swing.GroupLayout.PREFERRED_SIZE))
            .addGroup(imagePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(imagePanelLayout.createSequentialGroup()
                    .addComponent(image, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGap(29, 29, 29)))
        );

        controls.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(102, 102, 102)));

        jLabel2.setText("C:");

        jLabel3.setText("k:");

        jLabel4.setText("Speed:");

        input_constant.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                input_constantActionPerformed(evt);
            }
        });

        input_k.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                input_kActionPerformed(evt);
            }
        });

        jLabel5.setText("Settings:");

        buttonOK.setText("OK");
        buttonOK.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonOKActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout controlsLayout = new javax.swing.GroupLayout(controls);
        controls.setLayout(controlsLayout);
        controlsLayout.setHorizontalGroup(
            controlsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(controlsLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(controlsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel5, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(controlsLayout.createSequentialGroup()
                        .addGroup(controlsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel2)
                            .addComponent(jLabel3)
                            .addComponent(jLabel4))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addGroup(controlsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(input_speed, javax.swing.GroupLayout.DEFAULT_SIZE, 70, Short.MAX_VALUE)
                            .addComponent(input_k)
                            .addComponent(input_constant))))
                .addContainerGap())
            .addGroup(controlsLayout.createSequentialGroup()
                .addGap(45, 45, 45)
                .addComponent(buttonOK, javax.swing.GroupLayout.PREFERRED_SIZE, 70, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(47, Short.MAX_VALUE))
        );
        controlsLayout.setVerticalGroup(
            controlsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(controlsLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel5)
                .addGap(24, 24, 24)
                .addGroup(controlsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(input_constant, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(controlsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(input_k, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(controlsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel4)
                    .addComponent(input_speed, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(buttonOK)
                .addContainerGap(174, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addComponent(imagePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(controls, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                .addGap(0, 0, 0)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(controls, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addGap(30, 30, 30))
                    .addComponent(imagePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void input_constantActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_input_constantActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_input_constantActionPerformed

    private void input_kActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_input_kActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_input_kActionPerformed

    private void buttonOKActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonOKActionPerformed
        
        runAnimation();
    }//GEN-LAST:event_buttonOKActionPerformed

    private void imageComponentResized(java.awt.event.ComponentEvent evt) {//GEN-FIRST:event_imageComponentResized
        try{
            image.colorImage(image.getHeight(), image.getWidth(), u);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }//GEN-LAST:event_imageComponentResized

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(wave.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(wave.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(wave.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(wave.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new wave().setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton buttonOK;
    private javax.swing.JPanel controls;
    private mm.Image image;
    private javax.swing.JPanel imagePanel;
    private javax.swing.JTextField input_constant;
    private javax.swing.JTextField input_k;
    private javax.swing.JTextField input_speed;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JPanel jPanel1;
    // End of variables declaration//GEN-END:variables
}
