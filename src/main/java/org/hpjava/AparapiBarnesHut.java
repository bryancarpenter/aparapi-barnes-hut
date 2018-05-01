
package org.hpjava;

import java.awt.Color ;
import java.awt.Dimension ;
import java.awt.Graphics ;
import javax.swing.JFrame ;
import javax.swing.JPanel ;

import java.util.Random ;

import com.aparapi.Range ;
import com.aparapi.internal.kernel.KernelManager ;

public class AparapiBarnesHut {

    // Size of simulation

    //final static int N = 10000 ;  // Number of "stars"
    final static int N = 250000 ;  // Number of "stars"
    final static float BOX_WIDTH = 100.0F ;

    // Initial state

    final static float RADIUS = 20.0F ;  // of randomly populated sphere

    //final static float ANGULAR_VELOCITY = 0.4F ;
    //final static float ANGULAR_VELOCITY = 3F ;
    final static float ANGULAR_VELOCITY = 1.5F ;
           // controls total angular momentum (tend to increase this
           // as N increases, to keep "galaxy" stable).


    // Simulation

    final static float DT = 0.0005F ;  // Time step
           // (tend to decrease this as N increases, to maintain accuracy)


    // Display

    final static int WINDOW_SIZE = 1000 ;
    final static int DELAY = 0 ;
    final static int OUTPUT_FREQ = 1 ;


    // Star positions
    static float [] x = new float [N] ;
    static float [] y = new float [N] ;
    static float [] z = new float [N] ;

    // Star velocities
    static float [] vx = new float [N] ;
    static float [] vy = new float [N] ;
    static float [] vz = new float [N] ;

    // Star accelerations
    static float [] ax = new float [N] ;
    static float [] ay = new float [N] ;
    static float [] az = new float [N] ;

    // Barnes Hut tree
    static Node tree ;
    
    static Display display = new Display() ;
    
    public static void main(String args []) throws Exception {

        // Define initial state of stars

        Random rand = new Random(1234) ;

        // Randomly choose plane for net angular velocity

        double nx = 2 * rand.nextDouble() - 1 ;
        double ny = 2 * rand.nextDouble() - 1 ;
        double nz = 2 * rand.nextDouble() - 1 ;
        double norm = 1.0 / Math.sqrt(nx * nx + ny * ny + nz * nz) ;
        nx *= norm ;
        ny *= norm ;
        nz *= norm ;
       

        // ... or just rotate in x, y plane
        //double nx = 0, ny = 0, nz = 1.0 ;

        // ... or just rotate in x, z plane
        //double nx = 0, ny = 1.0, nz = 0 ;

        // Initial state is ball of randomly placed stars, rotating
        // according to angular velocity chosen above.

        for(int i = 0 ; i < N ; i++) {

            // Place star randomly in sphere of specified radius
            double rx, ry, rz, r ;
            do {
                rx = (2 * rand.nextDouble() - 1) * RADIUS ;
                ry = (2 * rand.nextDouble() - 1) * RADIUS ;
                rz = (2 * rand.nextDouble() - 1) * RADIUS ;
                r = Math.sqrt(rx * rx + ry * ry + rz * rz) ;
            } while(r > RADIUS) ;

            x [i] = (float) (0.5 * BOX_WIDTH + rx) ;
            y [i] = (float) (0.5 * BOX_WIDTH + ry) ;
            z [i] = (float) (0.5 * BOX_WIDTH + rz) ;

            vx [i] = (float) (ANGULAR_VELOCITY * (ny * rz - nz * ry)) ; 
            vy [i] = (float) (ANGULAR_VELOCITY * (nz * rx - nx * rz)) ; 
            vz [i] = (float) (ANGULAR_VELOCITY * (nx * ry - ny * rx)) ; 
        }
         
        int iter = 0 ;
        while(true) {
            double dtOver2 = 0.5 * DT;
            double dtSquaredOver2 = 0.5 * DT * DT;  

            if(iter % OUTPUT_FREQ == 0) {
                System.out.println("iter = " + iter + ", time = " + iter * DT) ;
                display.repaint() ;
            }

            // Verlet integration:
            // http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet


            for (int i = 0; i < N; i++) {
                // update position
                // mod implements periodic box
                x[i] = mod((float) (x [i] + (vx[i] * DT) +
                                   (ax[i] * dtSquaredOver2)), BOX_WIDTH);
                y[i] = mod((float) (y [i] + (vy[i] * DT) +
                                   (ay[i] * dtSquaredOver2)), BOX_WIDTH);
                z[i] = mod((float) (z [i] + (vz[i] * DT) +
                                   (az[i] * dtSquaredOver2)), BOX_WIDTH);
                // update velocity halfway
                vx[i] += (ax[i] * dtOver2);
                vy[i] += (ay[i] * dtOver2);
                vz[i] += (az[i] * dtOver2);
            }    

            computeAccelerations();

            for (int i = 0; i < N; i++) {
                // finish updating velocity with new acceleration
                vx[i] += (ax[i] * dtOver2);
                vy[i] += (ay[i] * dtOver2);
                vz[i] += (az[i] * dtOver2);
            }

            iter++ ;
            //break ; // debug
        }       


    }

    static boolean reported = false ;

    // Compute accelerations of all stars from current positions:
    static void computeAccelerations() {
        
       
        // Build the BH tree

        long startTreeTime = System.currentTimeMillis();

        Node.numNodes = 0 ;

        tree = new Node(BOX_WIDTH / 2, BOX_WIDTH / 2, BOX_WIDTH / 2,
                        BOX_WIDTH) ;
        for (int i = 0; i < N; i++) {
            tree.addParticle(x [i], y [i], z [i]);
        }

        long endTreeTime = System.currentTimeMillis();
        System.out.println("time to build Tree = " +
                           (endTreeTime - startTreeTime) + " milliseconds"); 
        
        // Prcomputations on BH tree - also allocate nodes of "flattened"
        // tree in Java arrays.

        KernelTree kernel = new KernelTree(x, y, z, ax, ay, az, Node.numNodes) ;

        long startPreComputeTime = System.currentTimeMillis();
        
        tree.preComputeAndAllocateFlat(kernel) ;

        long endPreComputeTime = System.currentTimeMillis();
        System.out.println("time to precompute Tree = " +
                           (endPreComputeTime - startPreComputeTime) + " milliseconds"); 
        
        // Finish initialization of flattened tree nodes (set pointers)

        // tree in Java arrays.
        long startFlattenTime = System.currentTimeMillis();

        tree.flatten(KernelTree.NULL, KernelTree.NULL, kernel) ;

        long endFlattenTime = System.currentTimeMillis();
        System.out.println("time to flatten Tree = " +
                           (endFlattenTime - startFlattenTime) + " milliseconds"); 

        System.out.println("Number of nodes = " + Node.numNodes);
       
        // Interaction forces (gravity)
        // This is where the program spends most of its time.

        long startForceTime = System.currentTimeMillis();

        kernel.execute(Range.create(N)) ;  // Invoke code on GPU
        
        long endForceTime = System.currentTimeMillis();
        System.out.println("time to calculate forces = " +
                           (endForceTime - startForceTime) + " milliseconds");

        // Report on execution mode...
        if(!reported) {
            StringBuilder builder = new StringBuilder() ;
            KernelManager.instance().reportDeviceUsage(builder, true) ;
            System.out.println(builder) ;
            reported = true ;
        }

        kernel.dispose() ;
    }
    
    static class Display extends JPanel {

        static final double SCALE = WINDOW_SIZE / BOX_WIDTH ;

        Display() {

            setPreferredSize(new Dimension(WINDOW_SIZE, WINDOW_SIZE)) ;

            JFrame frame = new JFrame("MD");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {
            g.setColor(Color.BLACK) ;
            g.fillRect(0, 0, WINDOW_SIZE, WINDOW_SIZE) ;
            g.setColor(Color.WHITE) ;
            for(int i = 0 ; i < N ; i++) {
                int gx = (int) (SCALE * x [i]) ;
                int gy = (int) (SCALE * y [i]) ;
                if(0 <= gx && gx < WINDOW_SIZE && 0 < gy && gy < WINDOW_SIZE) { 
                    g.fillRect(gx, gy, 1, 1) ;
                }
            } 
        }
    }
    
    static float mod(float x, float box) {
        float reduced = x - ((int) (x / box) * box) ;
        return reduced >= 0 ? reduced : reduced + box ;
    }


    static class Node {

        // Node of Barnes-Hut tree

        final static float OPENING_ANGLE = 1.0F ;

        float xMid, yMid, zMid ;
        float size ;
        int nParticles ;
        float xCent, yCent, zCent ;  // centre of mass
        Node [] children ;

        float threshold ;

        int flatNode ;  // index of corresponding elements in flattened tree

        static int numNodes ;

        Node(float xMid, float yMid, float zMid, float size) {
            this.xMid = xMid ;
            this.yMid = yMid ;
            this.zMid = zMid ;
            this.size = size ;

            numNodes++ ;
        }

        void addParticle(float x, float y, float z) {
            /* In single precision following test sometimes fails through rounding erros
            float sizeBy2 = size / 2 ;
            if(x < xMid - sizeBy2 || x > xMid + sizeBy2 ||
               y < yMid - sizeBy2 || y > yMid + sizeBy2 ||
               z < zMid - sizeBy2 || z > zMid + sizeBy2) {
                System.out.println("x = " + x + ", y = " + y + ", z = " + z) ;  // debug
                throw new IllegalArgumentException("particle position outside " +
                                                   "bounding box of Node") ;
            }
            */
            if(nParticles == 0) {
                xCent = x ;
                yCent = y ;
                zCent = z ;
                nParticles = 1 ;
                return ;
            } 
            if(nParticles == 1) {
                children = new Node [8] ;
                addParticleToChild(xCent, yCent, zCent) ;  
            }
            addParticleToChild(x, y, z) ;  
            nParticles++ ;
        }

        void addParticleToChild(float x, float y, float z) {

            int childIdx = ((x < xMid) ? 0 : 4) + ((y < yMid) ? 0 : 2) +
                           ((z < zMid) ? 0 : 1) ;

            Node child = children [childIdx] ;
            if(child == null) {
                float sizeBy4 = size / 4 ;
                child = new Node((x < xMid) ? xMid - sizeBy4 : xMid + sizeBy4,
                                 (y < yMid) ? yMid - sizeBy4 : yMid + sizeBy4,
                                 (z < zMid) ? zMid - sizeBy4 : zMid + sizeBy4,
                                 size / 2) ;
                children [childIdx] = child ;
            }
            child.addParticle(x, y, z) ;
        }


        void preComputeAndAllocateFlat(KernelTree kernel) {

            // Precompute Centre of Mass of this node (where non-leaf node)
            // and opening threshold for force calculation.

            // In this implementation, also allocate nodes of flattened tree.

            if(children != null) {
                float xSum = 0, ySum = 0, zSum = 0 ;

                for(int i = 0 ; i < 8 ; i++) {
                    Node child = children [i] ;
                    if(child != null) {
                        child.preComputeAndAllocateFlat(kernel) ;
                        int nChild = child.nParticles ;
                        xSum += nChild * child.xCent ;
                        ySum += nChild * child.yCent ;
                        zSum += nChild * child.zCent ;
                    }
                }
                xCent = xSum / nParticles ;
                yCent = ySum / nParticles ;
                zCent = zSum / nParticles ;
            }

            float delta = distance(xCent, yCent, zCent) ;
            threshold = size / OPENING_ANGLE + delta ;

            flatNode = kernel.allocateNode(xMid, yMid, zMid, nParticles,
                                           xCent, yCent, zCent, threshold) ;
                                          
        }

        void flatten(int parent, int next, KernelTree kernel) {

            // Set up pointers to related nodes in flattened node:

            //   parent = parent node
            //   firstChild = head of linked list of children of this node
            //   next = next sibling in parent's list of children

            int firstChild = KernelTree.NULL ;
            if (children != null) {
                Node prev = null ;
                for(int i = 0 ; i < 8 ; i++) {
                    Node child = children [i] ;
                    if(child != null) {
                        if(firstChild == KernelTree.NULL) 
                            firstChild = child.flatNode ;
                        if(prev != null)
                            prev.flatten(flatNode, child.flatNode, kernel) ;
                        prev = child ;
                    }
                }
                if(prev != null)
                    prev.flatten(flatNode, KernelTree.NULL, kernel) ;
            }

            kernel.setNeighbours(flatNode, parent, firstChild, next) ;
        }

        float distance(float x, float y, float z) {

            // Distance from mid-point of this node.

            // Implements cyclic boundaries in spatial cube

            float dx, dy, dz;  // separations in x and y directions
            float dx2, dy2, dz2, rSquared ;
            dx = x - xMid ;
            if(dx > BOX_WIDTH / 2) dx -= BOX_WIDTH ;
            if(dx < -BOX_WIDTH / 2) dx += BOX_WIDTH ;
            dy = y - yMid ;
            if(dy > BOX_WIDTH / 2) dy -= BOX_WIDTH ;
            if(dy < -BOX_WIDTH / 2) dy += BOX_WIDTH ;
            dz = z - zMid ;
            if(dz > BOX_WIDTH / 2) dz -= BOX_WIDTH ;
            if(dz < -BOX_WIDTH / 2) dz += BOX_WIDTH ;
            dx2 = dx * dx;
            dy2 = dy * dy;
            dz2 = dz * dz;
            rSquared = dx2 + dy2 + dz2 ;
            return (float) Math.sqrt(rSquared) ;
        }
    }
}

