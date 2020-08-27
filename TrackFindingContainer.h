#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <set>

const int NGRIDBINX = 10; //number of x grid for the container
const int NGRIDBINY = 10; //number of y grid for the container
vector<int> search_mod;   //id for modules that has search region
vector<double> search_x;  //x coordinate for the search region
vector<double> search_y;  //y coordinate for the search region
double search_size = 600.; //mm

//addtional cuts for the new track finder algorithm
const double TrackProjPosCut[2] = {1e9, 1e9};
const double FineTrackProjPosCut[2] = {1e9, 1e9};
const double FineTrackMaxSlope[2] = {1e9, 1e9};

struct GridHit{
    int layer;
    int module;
    int index;
    bool used;
    double x;
    double y;
    double z;
    double resox;
    double resoy;
    GridHit(int l, int m, int a, bool b, double x1, double y1, double z1, double ex, double ey) 
    : layer(l), module(m), index(a), used(b), x(x1), y(y1), z(z1), resox(ex), resoy (ey) {
    }
};

//each layer has one GridHitContainer that contains hits from all
//the modules of this layer

//TODO: should really write it as a class if we decide to use something
//similar
struct GridHitContainer{
    int layer;
    std::map<int, double> centerx;
    std::map<int, double> centery;
    std::map<int, double> sizex;
    std::map<int, double> sizey;
    std::map<int, double> binwx;
    std::map<int, double> binwy;
    std::map<int, double> startx;
    std::map<int, double> starty;
    
    //dimension of the map should be equal to the number of modules this
    //layer has. 
    std::map<int, std::vector<std::vector<std::vector<GridHit*>>>>  gridArray;
    std::map<int, std::vector<GridHit*>> hitArray;
    
    GridHitContainer(int l) {layer = l;}
    
    void AddModule(int modID){
        if (hitArray.find(modID) != hitArray.end() &&
            gridArray.find(modID) != gridArray.end()){
            std::cout<<"module "<<modID<<" already exist"<<endl;
            return;
        }
    
        std::vector<std::vector<std::vector<GridHit*>>> a;
        std::vector<GridHit*> b;
        b.reserve(200);
        a.resize(NGRIDBINX);
        for (unsigned int i=0; i<a.size(); i++){
            a[i].resize(NGRIDBINY);
            for (unsigned int j=0; j<a[i].size(); j++){
                a[i][j].reserve(100);
            }
        }
        gridArray[modID] = a;
        hitArray[modID]  = b;
    }
    
    void SetGeoInfo(int mod, double cx, double cy, double sx, double sy){
        centerx[mod] = cx;
        centery[mod] = cy;
        sizex[mod]   = sx;
        sizey[mod]   = sy;
        startx[mod]  = cx - sx/2.;
        starty[mod]  = cy - sy/2.;
        binwx[mod]   = sx / NGRIDBINX;
        binwy[mod]   = sy / NGRIDBINY;
    }
    
    void Clear(){
        for (map<int, std::vector<std::vector<std::vector<GridHit*>>>>::iterator it = gridArray.begin(); it != gridArray.end(); ++it){
            int modID = it->first;
            for (unsigned int i=0; i<gridArray[modID].size(); i++){
                for (unsigned int j=0; j<gridArray[modID][i].size(); j++){
                    gridArray[modID][i][j].clear();
                }
            }
            for (unsigned int i=0; i<hitArray[modID].size(); i++) delete hitArray[modID][i];
            hitArray[modID].clear();
        }
        centerx.clear();
        centery.clear();
        sizex.clear();
        sizey.clear();
        binwx.clear();
        binwy.clear();
        startx.clear();
        starty.clear();
    }
    
    void AddHit(int lay, int mod, int idx, double x, double y, double z, double resox, double resoy){
        int i = (x - startx[mod]) / binwx[mod];
        int j = (y - starty[mod]) / binwy[mod];
        if (i < NGRIDBINX && i >= 0 && j < NGRIDBINY && j >= 0){
            GridHit* thisHit = new GridHit(lay, mod, idx, false, x, y, z, resox, resoy);
            hitArray[mod].push_back(thisHit);
            gridArray[mod][i][j].push_back(thisHit);
        }
    }
    
    //find all the rectangular boxes that have partial overlap with the search region (NOTE: only work if
    //the two edges of the rectangle are parallel to the x-y axes)
    std::vector<vector<GridHit*>*> GetHitInRange (double x, double y, double r, int mod, int& nhits){
        std::vector<vector<GridHit*>*> a;
        a.clear();
        nhits = 0;
        
        int lowx  = ((x-r) - startx[mod]) / binwx[mod];
        int highx = ((x+r) - startx[mod]) / binwx[mod] + 0.5;
        int lowy  = ((y-r) - starty[mod]) / binwy[mod];
        int highy = ((y+r) - starty[mod]) / binwy[mod] + 0.5;

        //completely out of the range 
        if (lowx >= NGRIDBINX || highx < 0) return a;
        else if (lowy >= NGRIDBINY || highy < 0) return a;
        
        //limit the range if part of the circle is outside
        if (lowx < 0) lowx = 0;
        if (highx >= NGRIDBINX) highx = NGRIDBINX - 1;
        
        if (lowy < 0) lowy = 0;
        if (highy >= NGRIDBINY) highy = NGRIDBINY - 1;
        
        for (int i = lowx; i<= highx; i++){
            for (int j = lowy; j<= highy; j++){
                // Find the closest point to the circle within the rectangle
                double minx = startx[mod] + i*binwx[mod];
                double maxx = startx[mod] + (i+1)*binwx[mod];
                double miny = starty[mod] + j*binwy[mod];
                double maxy = starty[mod] + (j+1)*binwy[mod];
                
                double closestX = (x < minx) ? minx : (maxx < x) ? maxx : x;
                double closestY = (y < miny) ? miny : (maxy < y) ? maxy : y;

                // Calculate the distance between the circle's center and this closest point
                double dX = x - closestX;
                double dY = y - closestY;

                // If the distance is less than the circle's radius, an intersection occurs
                float d2 = (dX * dX) + (dY * dY);
                if (d2 < (r*r)) {
                    a.push_back(&gridArray[mod][i][j]);
                    nhits += gridArray[mod][i][j].size();
                }
            }
        }
        return a;
    }
    
    bool HasModule(int mod){
        std::map<int, std::vector<std::vector<std::vector<GridHit*>>>>::iterator it = gridArray.find(mod);
        if (it == gridArray.end()) return false;
        else return true;
    }
};

struct SBSReconTrack{
    bool status;
    std::vector<GridHit*> hits;
    double chi2ndf;
    double trackPara[4]; //x intersect, x slope, y intersect, y slope
    int nmissinghits;
    
    double sum[6];
    
    SBSReconTrack(int n) : nmissinghits (n) {
        status = true;
        hits.clear();
        chi2ndf = 0.;
        for (int i=0; i<4; i++) trackPara[i] = 0.;
        for (int i=0; i<6; i++) sum[i] = 0.;
    }
    
    //used when copying tracks
    SBSReconTrack(SBSReconTrack* a){
        status = a->status;
        chi2ndf= a->chi2ndf;
        hits   = a->hits;
        nmissinghits = a->nmissinghits;
        for (int i=0; i<4; i++) trackPara[i] = a->trackPara[i];
        for (int i=0; i<6; i++) sum[i] = a->sum[i];
    }
    
    void SetTrackPara(double a, double b, double c, double d){
        trackPara[0] = a;
        trackPara[1] = b;
        trackPara[2] = c;
        trackPara[3] = d;
    }
    
    void AddAndFilterWithLLS(GridHit* a){
        //update track parameters using simple linear least square (LLS) fit
        //right now using only the same method in the standalone code
        
        sum[0] += (a->x);
        sum[1] += (a->y);
        sum[2] += (a->z);
        sum[3] += (a->x)*(a->z);
        sum[4] += (a->y)*(a->z);
        sum[5] += pow((a->z),2);
        
        if (hits.size() == 1){
            //adding the second hit
            //now we can have inital estimate for track parameters
            trackPara[1] = (a->x - hits[0]->x) / (a->z - hits[0]->z);
            trackPara[3] = (a->y - hits[0]->y) / (a->z - hits[0]->z);
            
            trackPara[0] = hits[0]->x - trackPara[1]*hits[0]->z;
            trackPara[2] = hits[0]->y - trackPara[3]*hits[0]->z;
        }
        else if (hits.size() > 1){
            //adding the third hit and beyond
            int nhits = hits.size() + 1;
            double denom = (sum[5]*nhits - pow(sum[2],2));
            
            trackPara[0] = (sum[5]*sum[0] - sum[2]*sum[3])/denom;
            trackPara[1] = (nhits*sum[3] - sum[0]*sum[2])/denom;
            trackPara[2] = (sum[5]*sum[1] - sum[2]*sum[4])/denom;
		    trackPara[3] = (nhits*sum[4] - sum[1]*sum[2])/denom;
            
        }
        
        hits.push_back(a);
        
    }
    
    void CalChi2NDF(){
        int ndf = 2*hits.size()-4;
        double chi2 = 0.;
        for (unsigned int i=0; i<hits.size(); i++){
            double prox = trackPara[0] + trackPara[1]*(hits[i])->z;
            double proy = trackPara[2] + trackPara[3]*(hits[i])->z;
            chi2 += pow((prox - (hits[i])->x)/(hits[i])->resox, 2) 
                  + pow((proy - (hits[i])->y)/(hits[i])->resoy, 2);
        }
        chi2ndf = chi2 / ndf;
    }
    
    void GetProjection(double tarz, double& projx, double& projy){
        assert(hits.size() >= 2);
        projx = trackPara[0] + trackPara[1]*tarz;
        projy = trackPara[2] + trackPara[3]*tarz;
    }
    
    void CheckOutlier(double testr){
        //exclude each hit and fit the rest to check if the excluded hit
        //is compatible with the rest
        if (hits.size() <= 2) return;
        int index = -1;
        for (unsigned int i=0; i<hits.size(); i++){
            double csum[6] = {0};
            for (unsigned int j=0; j<hits.size(); j++){
                if (j == i) continue;
                csum[0] += (hits[j]->x);
                csum[1] += (hits[j]->y);
                csum[2] += (hits[j]->z);
                csum[3] += (hits[j]->x)*(hits[j]->z);
                csum[4] += (hits[j]->y)*(hits[j]->z);
                csum[5] += pow((hits[j]->z),2);
            }
            int nhits = hits.size() - 1; //excluded one hit
            double denom = (csum[5]*nhits - pow(csum[2],2));
            
            double tr_x  = (csum[5]*csum[0] - csum[2]*csum[3])/denom;
            double tr_tx = (nhits*csum[3] - csum[0]*csum[2])/denom;
            double tr_y  = (csum[5]*csum[1] - csum[2]*csum[4])/denom;
		    double tr_ty = (nhits*csum[4] - csum[1]*csum[2])/denom;
		    
		    double projx = tr_x + tr_tx * hits[i]->z;
		    double projy = tr_y + tr_ty * hits[i]->z;
		    
		    double dist  = sqrt(pow(projx - hits[i]->x, 2) 
		                      + pow(projy - hits[i]->y, 2));
		    if (dist > testr) index = i;
        }
        
        if (index >= 0) {
            hits.erase(hits.begin() + index);
            nmissinghits++;
        }
    }
    
    bool ForwardSearchTest(double testr){
        if (hits.size() <= 2) return false;
        //check also forward searching
        for (unsigned int i = 2; i<hits.size(); i++){
            double csum[6] = {0};
            for (unsigned int j=0; j<i; j++){
                csum[0] += (hits[j]->x);
                csum[1] += (hits[j]->y);
                csum[2] += (hits[j]->z);
                csum[3] += (hits[j]->x)*(hits[j]->z);
                csum[4] += (hits[j]->y)*(hits[j]->z);
                csum[5] += pow((hits[j]->z),2);
            }
            int nhits = i; 
            double denom = (csum[5]*nhits - pow(csum[2],2));
            
            double tr_x  = (csum[5]*csum[0] - csum[2]*csum[3])/denom;
            double tr_tx = (nhits*csum[3] - csum[0]*csum[2])/denom;
            double tr_y  = (csum[5]*csum[1] - csum[2]*csum[4])/denom;
		    double tr_ty = (nhits*csum[4] - csum[1]*csum[2])/denom;
		    
		    double projx = tr_x + tr_tx * hits[i]->z;
		    double projy = tr_y + tr_ty * hits[i]->z;
		    
		    double dist  = sqrt(pow(projx - hits[i]->x, 2) 
		                      + pow(projy - hits[i]->y, 2));
		                      
		    if (dist > testr) return true;
        }
        
        return false;
    }
        
};

//used for sorting recon tracks using c++ sort algorithm 
//so the first track should have the maximum number of hits,
// if there are multiple hits with the same amount the recon hits
//then the one with smaller chi2/ndf rank first
bool SortTracks(SBSReconTrack const& a, SBSReconTrack const& b){

    if (a.hits.size() > b.hits.size()) return true;
    else if (a.hits.size() == b.hits.size()){
        if (a.chi2ndf < b.chi2ndf) return true;
        else return false;
    }
    else return false;
    
}

bool SortTracksStatus(SBSReconTrack const& a, SBSReconTrack const& b){

    if (a.status > b.status) return true;
    else return false;
    
}

bool SortHits(GridHit* a, GridHit* b){
    return a->layer < b->layer;
}

