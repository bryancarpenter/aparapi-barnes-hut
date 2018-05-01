
package org.hpjava;

import com.aparapi.Kernel;


/*

To pass the tree through to kernel running on graphics device, BH
tree has to be "flattened" into array of Java primitive types.

Each original node is represented by NODEDSIZE consective elements
in a big float array, and NODEISIZE consecutive elements in a big int array,
with defined offsets for fields.

*/

public class KernelTree extends Kernel {


    final static float BOX_WIDTH = AparapiBarnesHut.BOX_WIDTH ;

    // Node structure constants
    
    // offsets of float fields
    
    final static int XMID = 0 ;
    final static int YMID = 1 ;
    final static int ZMID = 2 ;
    
    final static int XCENT = 3 ;
    final static int YCENT = 4 ;
    final static int ZCENT = 5 ;
    
    final static int THRESHOLD = 6 ;
    
    final static int NODEDSIZE = 7 ;
    
    // offsets of int fields
    
    final static int PARENT = 0 ;
    final static int FIRSTCHILD = 1 ;
    final static int NEXT = 2 ;
    
    final static int NPARTICLES = 3 ;
    
    final static int NODEISIZE = 4 ;
    
    
    final static int NULL = 0 ;        

    final static int TREE_ROOT = 1 ;  // assumed first allocated

        
    // Star positions
    final float [] x ;
    final float [] y ;
    final float [] z ;

    // Star accelerations
    final float [] ax ;
    final float [] ay ;
    final float [] az ;

    final float [] nodesD ;
    final int [] nodesI ;
    
    int nodeTop ;

    /* Constructor allocates big arrays.  By the time this is called
     * the size of the BH tree is known.
     */
    KernelTree(float [] x, float [] y, float [] z,
               float [] ax, float [] ay, float [] az, int numNodes) {
  
        this.x = x ;
        this.y = y ;
        this.z = z ;

        this.ax = ax ;
        this.ay = ay ;
        this.az = az ;

        int num = numNodes + 1 ;  // reserve 0 to represent NULL
        nodesD = new float [NODEDSIZE * num] ;
        nodesI = new int [NODEDSIZE * num] ;

        nodeTop = 1 ;  // 0 reserved
    }

    int allocateNode(float xMid, float yMid, float zMid, int nParticles,
                     float xCent, float yCent, float zCent,
                     float threshold) {

        /*
         * Allocate space for node and define most fields.
         */

        int node = nodeTop++ ;

        int nodesDptr = NODEDSIZE * node ;
        int nodesIptr = NODEISIZE * node ;

        nodesD [nodesDptr + XMID] = xMid ;
        nodesD [nodesDptr + YMID] = yMid ;
        nodesD [nodesDptr + ZMID] = zMid ;

        nodesI [nodesIptr + NPARTICLES] = nParticles ;

        nodesD [nodesDptr + XCENT] = xCent ;
        nodesD [nodesDptr + YCENT] = yCent ;
        nodesD [nodesDptr + ZCENT] = zCent ;

        nodesD [nodesDptr + THRESHOLD] = threshold ;

        return node ;
    }

    void setNeighbours(int node, int parent, int firstChild, int next) {

        /*
         * Set pointers to related nodes
         */

        int nodesIptr = NODEISIZE * node ;

        nodesI [nodesIptr + PARENT] = parent ;
        nodesI [nodesIptr + FIRSTCHILD] = firstChild ;
        nodesI [nodesIptr + NEXT] = next;
    }


    /*
     *
     * For reference below: depth first traversal without recursion...
     * Here children of tree represented by a linked list (firstChild
     * and next fields).
     *
     * Node current = root
     * boolean firstVisit = true ;
     *
     * while(true) {
     *    if(firstVisit) {
     *        if(current.firstChild != null) {
     *            ... processing of current node before visiting children ...
     *            current = current.firstChild ;
     *            continue ;
     *        }
     *        else {
     *            ... processing of current node if a leaf node ...
     *        }
     *    }
     *    else {
     *        ... processing of current node after visiting children ...
     *    }
     *    if(current.next != null) {
     *        current = current.next ;  // next sibling node
     *        firstVisit = true ;
     *        continue ;
     *    }
     *    if(current.parent != null) {
     *        current = current.parent ;
     *        firstVisit = false ;
     *        continue ;
     *    }
     *    break ;
     * }
     */

    void calcForce(int id, float x, float y, float z, int tree) {

        /*
         * Has to be rewritten because GPU kernels don't generally
         * support recursion.  Also there are a few quirks of Aparapi,
         * in terms of what statements it can generate GPU code for.
         */
               
        int current = tree ;
        boolean firstVisit = true ;

        int nodesIptr ;       // Aparapi doesn't like declarations inside while
        boolean continuing ;  // Aparapi doesn't like continue statement
        boolean breaking = false ;  // Aparapi fussy about break statements
        while(!breaking) {
            continuing = false ;

            nodesIptr = NODEISIZE * current ;            
            if(firstVisit) {
                if(!calcForceNodeRule(id, x, y, z, current)) {
                //if(nodesI [nodesIptr + FIRSTCHILD] != NULL) {  // debug
                    // visit children
                    current = nodesI [nodesIptr + FIRSTCHILD] ;
                    continuing = true ;
                }
            }

            if(!continuing) {
                int next = nodesI [nodesIptr + NEXT] ;
                if(next != NULL) {
                    current = next ;
                    firstVisit = true ;
                    continuing = true ;
                }
            }

            if(!continuing) {
                int parent = nodesI [nodesIptr + PARENT] ;                
                if(parent != NULL) {
                    current = parent ;
                    firstVisit = false ;
                    continuing = true ; ;
                }
            }

            if(!continuing)
                breaking = true ;
        }            
    }
    
    boolean calcForceNodeRule(int id, float x, float y, float z,
                              int node) {
    
        // Return value of true means processing completed with this node.
        // Return value of false means children must be visited.
        
        int nodesDptr = NODEDSIZE * node ;
        int nodesIptr = NODEISIZE * node ;
        
/* Aparapi doesn't like new
        if(nodesI [nodesIptr + NPARTICLES] == 0)
            throw new RuntimeException("Node without any particles") ;
*/
        if(nodesI [nodesIptr + FIRSTCHILD] == NULL) {
            // leaf node
            if(x != nodesD [nodesDptr + XCENT] || 
               y != nodesD [nodesDptr + YCENT] ||
               z != nodesD [nodesDptr + ZCENT]) {
                forceLaw(id, x, y, z, node) ;
            }
            return true ;
        }
        else {                 
            float r = distance(x, y, z, node) ;
            if(r > nodesD [nodesDptr + THRESHOLD]) {
                forceLaw(id, x, y, z, node) ;
                return true ;
            }
            else {
                return false ;
            }
            //return false ;  // debug
        }
    }
    

    float distance(float x, float y, float z, int node) {

        // Distance from mid-point of this node.

        int nodesDptr = NODEDSIZE * node ;
        
        float dx, dy, dz;  // separations in x and y directions
        float dx2, dy2, dz2, rSquared ;
        dx = x - nodesD [nodesDptr + XMID] ;
        if(dx > BOX_WIDTH / 2) dx -= BOX_WIDTH ;
        if(dx < -BOX_WIDTH / 2) dx += BOX_WIDTH ;
        dy = y - nodesD [nodesDptr + YMID] ;
        if(dy > BOX_WIDTH / 2) dy -= BOX_WIDTH ;
        if(dy < -BOX_WIDTH / 2) dy += BOX_WIDTH ;
        dz = z - nodesD [nodesDptr + ZMID] ;
        if(dz > BOX_WIDTH / 2) dz -= BOX_WIDTH ;
        if(dz < -BOX_WIDTH / 2) dz += BOX_WIDTH ;
        dx2 = dx * dx;
        dy2 = dy * dy;
        dz2 = dz * dz;
        rSquared = dx2 + dy2 + dz2 ;
        //return (float) Math.sqrt(rSquared) ;
        return sqrt(rSquared) ;  // Aparapi implements sqrt natively
    }

    void forceLaw(int id, float x, float y, float z, int node) {

        // Force exerted by effective mass at CM of this node.
        
        int nodesDptr = NODEDSIZE * node ;
        int nodesIptr = NODEISIZE * node ;
        
        float dx, dy, dz;  // separations in x and y directions
        float dx2, dy2, dz2, rSquared, r, massRCubedInv;      

        // Vector version of inverse square law
        // This version assumes periodic box.
        dx = x - nodesD [nodesDptr + XCENT] ;
        if(dx > BOX_WIDTH / 2) dx -= BOX_WIDTH ;
        if(dx < -BOX_WIDTH / 2) dx += BOX_WIDTH ;
        dy = y - nodesD [nodesDptr + YCENT] ;
        if(dy > BOX_WIDTH / 2) dy -= BOX_WIDTH ;
        if(dy < -BOX_WIDTH / 2) dy += BOX_WIDTH ;
        dz = z - nodesD [nodesDptr + ZCENT] ;
        if(dz > BOX_WIDTH / 2) dz -= BOX_WIDTH ;
        if(dz < -BOX_WIDTH / 2) dz += BOX_WIDTH ;
        dx2 = dx * dx;
        dy2 = dy * dy;
        dz2 = dz * dz;
        rSquared = dx2 + dy2 + dz2 ;
        //r = (float) Math.sqrt(rSquared) ;
        r = sqrt(rSquared) ;  // Aparapi implements this natively
        massRCubedInv = nodesI [nodesIptr + NPARTICLES] / (rSquared * r) ;
        ax [id] -= massRCubedInv * dx ;
        ay [id] -= massRCubedInv * dy ;
        az [id] -= massRCubedInv * dz ;
    }

    public void run() {

        int gid = getGlobalId() ;

        ax [gid] = 0F ;
        ay [gid] = 0F ;
        az [gid] = 0F ;

        calcForce(gid, x [gid], y [gid], z [gid], TREE_ROOT) ;
    }
}

