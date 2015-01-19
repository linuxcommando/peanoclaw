#if defined(SWE) || defined(PEANOCLAW_FULLSWOF2D)

#include <algorithm>

#include <png.h>

#include "peanoclaw/native/scenarios/MekkaFlood.h"
#include "gsl/gsl_interp.h"

//#define BREAKINGDAMTEST

/*MekkaFlood_SWEKernelScenario::MekkaFlood_SWEKernelScenario(double _domainSize) :
    _domainSize(_domainSize),
    bathymetryHelper("smith_sandwell.nc", "X2MIN", "Y2MIN","ROSE")
{

}*/

peanoclaw::native::scenarios::MekkaFloodSWEScenario::MekkaFloodSWEScenario(
    const tarch::la::Vector<DIMENSIONS,int>&    subdivisionFactor,
    const tarch::la::Vector<DIMENSIONS,double>& minimalMeshWidth,
    const tarch::la::Vector<DIMENSIONS,double>& maximalMeshWidth,
    double                                      globalTimestepSize,
    double                                      endTime
) : _subdivisionFactor(subdivisionFactor),
    _minialMeshWidth(minimalMeshWidth),
    _maximalMeshWidth(maximalMeshWidth),
    _globalTimestepSize(globalTimestepSize),
    _endTime(endTime)
{

    dem.load("DEM_100cm.bin");
 
    scale = 2.0;

#if 0
    // open png with mekka_map
    memset(&mekka_map, 0, (sizeof mekka_map));
    mekka_map.version = PNG_IMAGE_VERSION;
 
    /* The first argument is the file to read: */
    if (png_image_begin_read_from_file(&mekka_map, "tex2d_cropped.png")) {
         /* Set the format in which to read the PNG file; this code chooses a
          * simple sRGB format with a non-associated alpha channel, adequate to
          * store most images.
          */
         //mekka_map.format = PNG_FORMAT_RGBA;
         mekka_map.format = PNG_FORMAT_GRAY; // we only need gray scale

         /* Now allocate enough memory to hold the image in this format; the
          * PNG_IMAGE_SIZE macro uses the information about the image (width,
          * height and format) stored in 'image'.
          */
         std::cout << "image buffer size: " << PNG_IMAGE_SIZE(mekka_map) << std::endl;
         mekka_map_data = static_cast<uint8_t*>(malloc(PNG_IMAGE_SIZE(mekka_map)));

         if (mekka_map_data != NULL &&
            png_image_finish_read(&mekka_map, NULL/*background*/, mekka_map_data,
               0/*row_stride*/, NULL/*colormap*/)) {
         
              std::cout << "got png image data: width=" << mekka_map.width << ", height=" << mekka_map.height << std::endl;
         }

    } else {

        /* Something went wrong reading or writing the image.  libpng stores a
         * textual message in the 'png_image' structure:
         */
         std::cerr << "pngtopng: error: " << mekka_map.message << std::endl;
    }
#endif

    double upper_right_0 = dem.upper_right(0);
    double upper_right_1 = dem.upper_right(1);

    double lower_left_0 = dem.lower_left(0);
    double lower_left_1 = dem.lower_left(1);

    // we have to divide because otherwise Peano gets stuck
    double x_size = (upper_right_0 - lower_left_0) / scale;
    double y_size = (upper_right_1 - lower_left_1) / scale;

    // TODO: make central scale parameter in MekkaFlood class
    // currently we have to change here, meshToCoordinates and initializePatch and computeMeshWidth
    assignList(_domainSize) = x_size, y_size;

    // TODO: aaarg Y U NO PLOT CORRECTLY! -> work around established
    _domainOffset = tarch::la::Vector<DIMENSIONS, double>(0);
    
}

peanoclaw::native::scenarios::MekkaFloodSWEScenario::~MekkaFloodSWEScenario() {
    // cleanup png data
#if 0
    if (mekka_map_data == NULL)
       png_image_free(&mekka_map);

    else
       free(mekka_map_data);
#endif
}

static double bilinear_interpolate(double x00, double y00,
                                   double x11, double y11,
                                   double x, double y,
                                   double f00, double f10, double f01, double f11) 
{
    // NOTE: dataset is organized as latitude / longitude

    double width = x11 - x00;
    double height = y11 - y00;

    // 1D interpolation along bottom
    double bottom_interpolated = ((x11 - x) * f00) / width + 
                                 ((x - x00) * f10) / width;

    // 1D interpolation along top
    double top_interpolated = ((x11 - x) * f01) / width + 
                              ((x - x00) * f11) / width;

    // complete the bilinear interpolation along latitude
    double data_value = ((y11 - y) * bottom_interpolated) / height + 
                        ((y - y00) * top_interpolated) / height;
 
    return data_value;
}

static double interpolation_error(peanoclaw::Patch& patch, int unknown) {
    tarch::la::Vector<DIMENSIONS, int> subcellIndex;
    const tarch::la::Vector<DIMENSIONS, double> patchSize = patch.getSize();
//    const tarch::la::Vector<DIMENSIONS, double> patchPosition = patch.getPosition();
    const tarch::la::Vector<DIMENSIONS, double> meshWidth = patch.getSubcellSize();
    const tarch::la::Vector<DIMENSIONS, int> subdivisionFactor = patch.getSubdivisionFactor();
 
    // interpolatation error inside patch
    double max_error = 0.0;
  
    subcellIndex(0) = 0;
    subcellIndex(1) = 0;
    double f00 = patch.getAccessor().getValueUNew(subcellIndex, 0);
 
    subcellIndex(0) = subdivisionFactor(0)-1;
    subcellIndex(1) = 0;
    double f10 = patch.getAccessor().getValueUNew(subcellIndex, 0);
 
    subcellIndex(0) = 0;
    subcellIndex(1) = subdivisionFactor(1)-1;
    double f01 = patch.getAccessor().getValueUNew(subcellIndex, 0);
 
    subcellIndex(0) = subdivisionFactor(0)-1;
    subcellIndex(1) = subdivisionFactor(1)-1;
    double f11 = patch.getAccessor().getValueUNew(subcellIndex, 0);

    for (int yi = 0; yi < subdivisionFactor(1); yi++) {
        for (int xi = 0; xi < subdivisionFactor(0); xi++) {
            subcellIndex(0) = xi;
            subcellIndex(1) = yi;
            double value = patch.getAccessor().getValueUNew(subcellIndex, 0);
 
            double interpolated_value = bilinear_interpolate(0.0, 0.0,
                                                             patchSize(0), patchSize(1),
                                                             xi*meshWidth(0), yi*meshWidth(0),
                                                             f00, f10, f01, f11);


            double error = std::abs(value - interpolated_value);

            max_error = std::max(max_error, error);
        }
    }

 
    // interpolatation error with ghostlayer
    double max_error_ghost = 0.0;
  
    subcellIndex(0) = -1;
    subcellIndex(1) = -1;
    f00 = patch.getAccessor().getValueUOld(subcellIndex, 0);
 
    subcellIndex(0) = subdivisionFactor(0);
    subcellIndex(1) = -1;
    f10 = patch.getAccessor().getValueUOld(subcellIndex, 0);
 
    subcellIndex(0) = -1;
    subcellIndex(1) = subdivisionFactor(1);
    f01 = patch.getAccessor().getValueUOld(subcellIndex, 0);
 
    subcellIndex(0) = subdivisionFactor(0);
    subcellIndex(1) = subdivisionFactor(1);
    f11 = patch.getAccessor().getValueUOld(subcellIndex, 0);

    for (int yi = -1; yi < subdivisionFactor(1)+1; yi++) {
        for (int xi = -1; xi < subdivisionFactor(0)+1; xi++) {
            subcellIndex(0) = xi;
            subcellIndex(1) = yi;
            double value = patch.getAccessor().getValueUOld(subcellIndex, 0);
 
            double interpolated_value = bilinear_interpolate(-meshWidth(0), -meshWidth(1),
                                                             patchSize(0)+meshWidth(0), patchSize(1)+meshWidth(1),
                                                             xi*meshWidth(0), yi*meshWidth(1),
                                                             f00, f10, f01, f11);


            double error = std::abs(value - interpolated_value);

            max_error_ghost = std::max(max_error, error);
        }
    }

    return std::max(max_error,max_error_ghost);
}


static double interpolation_error_coarse_fine_gradient(peanoclaw::Patch& patch, int unknown) {
    tarch::la::Vector<DIMENSIONS, int> subcellIndex;
    tarch::la::Vector<DIMENSIONS, int> coarseSubcellIndex;

//    const tarch::la::Vector<DIMENSIONS, double> patchSize = patch.getSize();
//    const tarch::la::Vector<DIMENSIONS, double> patchPosition = patch.getPosition();
    const tarch::la::Vector<DIMENSIONS, double> meshWidth = patch.getSubcellSize();
    const tarch::la::Vector<DIMENSIONS, int> subdivisionFactor = patch.getSubdivisionFactor();
 
    double max_error = 0.0;
 
    tarch::la::Vector<DIMENSIONS, double> meshPos;
    for (int yi = 1; yi < subdivisionFactor(1)-1; yi++) {
        for (int xi = 1; xi < subdivisionFactor(0)-1; xi++) {
            // coarse info ------------------------------
            coarseSubcellIndex(0) = std::floor(xi / 3.0) * 3;
            coarseSubcellIndex(1) = std::floor(yi / 3.0) * 3;
            tarch::la::Vector<DIMENSIONS, double> meshPos_00 = patch.getSubcellPosition(coarseSubcellIndex);
            double f00 = patch.getAccessor().getValueUNew(coarseSubcellIndex, unknown);
         
            coarseSubcellIndex(0) = std::ceil(xi / 3.0) * 3;
            coarseSubcellIndex(1) = std::floor(yi / 3.0) * 3;
            double f10 = patch.getAccessor().getValueUNew(coarseSubcellIndex, unknown);
         
            coarseSubcellIndex(0) = std::floor(xi / 3.0) * 3;
            coarseSubcellIndex(1) = std::ceil(yi / 3.0) * 3;
            double f01 = patch.getAccessor().getValueUNew(coarseSubcellIndex, unknown);
         
            coarseSubcellIndex(0) = std::ceil(xi / 3.0) * 3;
            coarseSubcellIndex(1) = std::ceil(yi / 3.0) * 3;
            tarch::la::Vector<DIMENSIONS, double> meshPos_11 = patch.getSubcellPosition(coarseSubcellIndex);
            double f11 = patch.getAccessor().getValueUNew(coarseSubcellIndex, unknown);

            // fine info -------------------------------
            // center
            subcellIndex(0) = xi;
            subcellIndex(1) = yi;
            meshPos = patch.getSubcellPosition(subcellIndex);
            double value_11 = patch.getAccessor().getValueUNew(subcellIndex, 0);
            double interpolated_value_11 = bilinear_interpolate(meshPos_00(0), meshPos_00(1),
                                                             meshPos_11(0), meshPos_11(1),
                                                             meshPos(0), meshPos(1),
                                                             f00, f10, f01, f11);
            double error_11 = value_11 - interpolated_value_11;
            
            // left
            subcellIndex(0) = xi-1;
            subcellIndex(1) = yi;
            meshPos = patch.getSubcellPosition(subcellIndex);
            double value_01 = patch.getAccessor().getValueUNew(subcellIndex, 0);
            double interpolated_value_01 = bilinear_interpolate(meshPos_00(0), meshPos_00(1),
                                                             meshPos_11(0), meshPos_11(1),
                                                             meshPos(0), meshPos(1),
                                                             f00, f10, f01, f11);
            double error_01 = value_01 - interpolated_value_01;

            // right
            subcellIndex(0) = xi+1;
            subcellIndex(1) = yi;
            meshPos = patch.getSubcellPosition(subcellIndex);
            double value_21 = patch.getAccessor().getValueUNew(subcellIndex, 0);
            double interpolated_value_21 = bilinear_interpolate(meshPos_00(0), meshPos_00(1),
                                                             meshPos_11(0), meshPos_11(1),
                                                             meshPos(0), meshPos(1),
                                                             f00, f10, f01, f11);
            double error_21 = value_21 - interpolated_value_21;

            // bottom
            subcellIndex(0) = xi;
            subcellIndex(1) = yi-1;
            meshPos = patch.getSubcellPosition(subcellIndex);
//            double value_10 = patch.getAccessor().getValueUNew(subcellIndex, 0);
//            double interpolated_value_10 = bilinear_interpolate(meshPos_00(0), meshPos_00(1),
//                                                             meshPos_11(0), meshPos_11(1),
//                                                             meshPos(0), meshPos(1),
//                                                             f00, f10, f01, f11);
////            double error_10 = value_10 - interpolated_value_10;
 
            // top
            subcellIndex(0) = xi;
            subcellIndex(1) = yi+1;
            meshPos = patch.getSubcellPosition(subcellIndex);
//            double value_12 = patch.getAccessor().getValueUNew(subcellIndex, 0);
//            double interpolated_value_12 = bilinear_interpolate(meshPos_00(0), meshPos_00(1),
//                                                             meshPos_11(0), meshPos_11(1),
//                                                             meshPos(0), meshPos(1),
//                                                             f00, f10, f01, f11);
//            double error_12 = value_12 - interpolated_value_12;


            double left_gradient = std::abs(error_11 - error_01);
            double right_gradient = std::abs(error_21 - error_11);

//            double bottom_gradient = std::abs(error_11 - error_10);
//            double top_gradient = std::abs(error_12 - error_11);


            double max_gradient_x = std::max(left_gradient, right_gradient) / meshWidth(0);
            double max_gradient_y = std::max(left_gradient, right_gradient) / meshWidth(1);

            max_error = std::max(max_error, max_gradient_x);
            max_error = std::max(max_error, max_gradient_y);
        }
    }

    return max_error;
}


static double interpolation_error_gradient(peanoclaw::Patch& patch, int unknown) {
    tarch::la::Vector<DIMENSIONS, int> subcellIndex;
//    const tarch::la::Vector<DIMENSIONS, double> patchSize = patch.getSize();
//    const tarch::la::Vector<DIMENSIONS, double> patchPosition = patch.getPosition();
    const tarch::la::Vector<DIMENSIONS, double> meshWidth = patch.getSubcellSize();
    const tarch::la::Vector<DIMENSIONS, int> subdivisionFactor = patch.getSubdivisionFactor();
 
    // interpolatation error inside patch
    double max_error = 0.0;
  
    subcellIndex(0) = 0;
    subcellIndex(1) = 0;
    tarch::la::Vector<DIMENSIONS, double> meshPos_00 = patch.getSubcellPosition(subcellIndex);
    double f00 = patch.getAccessor().getValueUNew(subcellIndex, unknown);
 
    subcellIndex(0) = subdivisionFactor(0)-1;
    subcellIndex(1) = 0;
    double f10 = patch.getAccessor().getValueUNew(subcellIndex, unknown);
 
    subcellIndex(0) = 0;
    subcellIndex(1) = subdivisionFactor(1)-1;
    double f01 = patch.getAccessor().getValueUNew(subcellIndex, unknown);
 
    subcellIndex(0) = subdivisionFactor(0)-1;
    subcellIndex(1) = subdivisionFactor(1)-1;
    tarch::la::Vector<DIMENSIONS, double> meshPos_11 = patch.getSubcellPosition(subcellIndex);
    double f11 = patch.getAccessor().getValueUNew(subcellIndex, unknown);
 
    tarch::la::Vector<DIMENSIONS, double> meshPos;
    for (int yi = 1; yi < subdivisionFactor(1)-1; yi++) {
        for (int xi = 1; xi < subdivisionFactor(0)-1; xi++) {
            // center
            subcellIndex(0) = xi;
            subcellIndex(1) = yi;
            meshPos = patch.getSubcellPosition(subcellIndex);
            double value_11 = patch.getAccessor().getValueUNew(subcellIndex, unknown);
            double interpolated_value_11 = bilinear_interpolate(meshPos_00(0), meshPos_00(1),
                                                             meshPos_11(0), meshPos_11(1),
                                                             meshPos(0), meshPos(1),
                                                             f00, f10, f01, f11);
            double error_11 = value_11 - interpolated_value_11;
            
            // left
            subcellIndex(0) = xi-1;
            subcellIndex(1) = yi;
            meshPos = patch.getSubcellPosition(subcellIndex);
            double value_01 = patch.getAccessor().getValueUNew(subcellIndex, unknown);
            double interpolated_value_01 = bilinear_interpolate(meshPos_00(0), meshPos_00(1),
                                                             meshPos_11(0), meshPos_11(1),
                                                             meshPos(0), meshPos(1),
                                                             f00, f10, f01, f11);
            double error_01 = value_01 - interpolated_value_01;

            // right
            subcellIndex(0) = xi+1;
            subcellIndex(1) = yi;
            meshPos = patch.getSubcellPosition(subcellIndex);
            double value_21 = patch.getAccessor().getValueUNew(subcellIndex, unknown);
            double interpolated_value_21 = bilinear_interpolate(meshPos_00(0), meshPos_00(1),
                                                             meshPos_11(0), meshPos_11(1),
                                                             meshPos(0), meshPos(1),
                                                             f00, f10, f01, f11);
            double error_21 = value_21 - interpolated_value_21;

            // bottom
            subcellIndex(0) = xi;
            subcellIndex(1) = yi-1;
            meshPos = patch.getSubcellPosition(subcellIndex);
//            double value_10 = patch.getAccessor().getValueUNew(subcellIndex, unknown);
//            double interpolated_value_10 = bilinear_interpolate(meshPos_00(0), meshPos_00(1),
//                                                             meshPos_11(0), meshPos_11(1),
//                                                             meshPos(0), meshPos(1),
//                                                             f00, f10, f01, f11);
//            double error_10 = value_10 - interpolated_value_10;
 
            // top
            subcellIndex(0) = xi;
            subcellIndex(1) = yi+1;
            meshPos = patch.getSubcellPosition(subcellIndex);
//            double value_12 = patch.getAccessor().getValueUNew(subcellIndex, 0);
//            double interpolated_value_12 = bilinear_interpolate(meshPos_00(0), meshPos_00(1),
//                                                             meshPos_11(0), meshPos_11(1),
//                                                             meshPos(0), meshPos(1),
//                                                             f00, f10, f01, f11);
//            double error_12 = value_12 - interpolated_value_12;


            double left_gradient = std::abs(error_11 - error_01);
            double right_gradient = std::abs(error_21 - error_11);

//            double bottom_gradient = std::abs(error_11 - error_10);
//            double top_gradient = std::abs(error_12 - error_11);


            double max_gradient_x = std::max(left_gradient, right_gradient) / meshWidth(0);
            double max_gradient_y = std::max(left_gradient, right_gradient) / meshWidth(1);

            max_error = std::max(max_error, max_gradient_x);
            max_error = std::max(max_error, max_gradient_y);
        }
    }

    return max_error;
}

void peanoclaw::native::scenarios::MekkaFloodSWEScenario::initializePatch(peanoclaw::Patch& patch) {
    // compute from mesh data
//    const tarch::la::Vector<DIMENSIONS, double> patchSize = patch.getSize();
//    const tarch::la::Vector<DIMENSIONS, double> patchPosition = patch.getPosition();
//    const tarch::la::Vector<DIMENSIONS, double> meshWidth = patch.getSubcellSize();
 
    const tarch::la::Vector<DIMENSIONS, double> tree_position = patch.getPosition();
    const tarch::la::Vector<DIMENSIONS, double> patchSize = patch.getSize();

    ssize_t domain_pixelspace_offset[2] = {round(tree_position[0]*scale),round(tree_position[1]*scale)};
  
    const tarch::la::Vector<DIMENSIONS, int> subdivisionFactor = patch.getSubdivisionFactor();
    ssize_t patch_cells[2];
    patch_cells[0] = subdivisionFactor[0];
    patch_cells[1] = subdivisionFactor[1];

    ssize_t patch_vertices[2];
    patch_vertices[0] = patch_cells[0]+1;
    patch_vertices[1] = patch_cells[1]+1;

    double local_pixelspace_size[2];
    local_pixelspace_size[0] = (patchSize[0]*scale) / (patch_vertices[0]-1);
    local_pixelspace_size[1] = (patchSize[1]*scale) / (patch_vertices[1]-1);


    // initialize new part only
    tarch::la::Vector<DIMENSIONS, int> subcellIndex;
    tarch::la::Vector<DIMENSIONS, double> meshPos;
    for (int yi = 0; yi < patch.getSubdivisionFactor()(1); yi++) {
        for (int xi = 0; xi < patch.getSubdivisionFactor()(0); xi++) {
            // evaluate cell centers
            ssize_t position[2] = {xi, yi};
 
            int pixelspace_x[2];
            int pixelspace_y[2];

            pixelspace_x[0] = domain_pixelspace_offset[0] + round(local_pixelspace_size[0] * position[0]);
            pixelspace_x[1] = domain_pixelspace_offset[0] + round(local_pixelspace_size[0] * (position[0]+1));
            pixelspace_y[0] = domain_pixelspace_offset[1] + round(local_pixelspace_size[1] * position[1]);
            pixelspace_y[1] = domain_pixelspace_offset[1] + round(local_pixelspace_size[1] * (position[1]+1));

            // cell corners
            float dem_height_x0y0 = 0.0;
            float dem_height_x1y0 = 0.0;
            float dem_height_x0y1 = 0.0;
            float dem_height_x1y1 = 0.0;

            if (pixelspace_x[0] >= 0 && pixelspace_x[1] < dem.dimension(0) 
             && pixelspace_y[0] >= 0 && pixelspace_y[1] <= dem.dimension(1)) 
            {
                dem_height_x0y0 = dem(pixelspace_x[0], pixelspace_y[0]);
                dem_height_x1y0 = dem(pixelspace_x[1], pixelspace_y[0]);
                dem_height_x0y1 = dem(pixelspace_x[0], pixelspace_y[1]);
                dem_height_x1y1 = dem(pixelspace_x[1], pixelspace_y[1]);
            }
 
#if !defined(NDEBUG)
            std::cout << "position: " << position[0] << " " << position[1] 
                      << " pixelspace x range: " << pixelspace_x[0] << " - " << pixelspace_x[1]-1
                      << " pixelspace y range: " << pixelspace_y[0] << " - " << pixelspace_y[1]-1
                      << " worldspace x range: " << pixelspace_x[0] * dem.scale(0) << " - " << (pixelspace_x[1]-1) * dem.scale(0)
                      << " worldspace y range: " << pixelspace_y[0] * dem.scale(1) << " - " << (pixelspace_y[1]-1) * dem.scale(1)
                      << " heights: " << dem_height_x0y0 << " " << dem_height_x1y0
                      << " " << dem_height_x0y1 << " " << dem_height_x1y1
                      << std::endl;
#endif

            double worldspace_x[2];
            double worldspace_y[2];

            worldspace_x[0] = pixelspace_x[0] * dem.scale(0);
            worldspace_x[1] = pixelspace_x[1] * dem.scale(0);
            worldspace_y[0] = pixelspace_y[0] * dem.scale(1);
            worldspace_y[1] = pixelspace_y[1] * dem.scale(1);

            double lower_heights[2] = {dem_height_x0y0, dem_height_x1y0};
            double upper_heights[2] = {dem_height_x0y1, dem_height_x1y1};

            // setup interpolation using cell corners
            gsl_interp *interp_lower = gsl_interp_alloc(gsl_interp_linear,2);
            gsl_interp_accel *allocp_lower = gsl_interp_accel_alloc();

            gsl_interp *interp_upper = gsl_interp_alloc(gsl_interp_linear,2);
            gsl_interp_accel *allocp_upper = gsl_interp_accel_alloc();

            gsl_interp_init(interp_lower, worldspace_x, lower_heights, 2);
            gsl_interp_init(interp_upper, worldspace_x, upper_heights, 2);

            gsl_interp *interp_temp = gsl_interp_alloc(gsl_interp_linear,2);
            gsl_interp_accel *allocp_temp = gsl_interp_accel_alloc();


            ssize_t dem_vertices[2];
            dem_vertices[0] = pixelspace_x[1] - pixelspace_x[0];
            dem_vertices[1] = pixelspace_y[1] - pixelspace_y[0];

            // interpolate on cell center going from lower left corner of cell
            int test_pixelspace_x = pixelspace_x[0];
            int test_pixelspace_y = pixelspace_y[0];

            double test_worldspace_x = (test_pixelspace_x + 0.5 * local_pixelspace_size[0]) * dem.scale(0);
            double test_worldspace_y = (test_pixelspace_y + 0.5 * local_pixelspace_size[1]) * dem.scale(1);

            double temp[2];
            temp[0] = gsl_interp_eval(interp_lower, worldspace_x, lower_heights, test_worldspace_x, allocp_lower);
            temp[1] = gsl_interp_eval(interp_upper, worldspace_x, upper_heights, test_worldspace_x, allocp_upper);

            gsl_interp_init(interp_temp, worldspace_y, temp, 2);
            double dem_interpolated = gsl_interp_eval(interp_temp, worldspace_y, temp, test_worldspace_y, allocp_temp);

            // assign set of values
            double h = 0.0;
            double u = 0.0;
            double v = 0.0;
            double mapvalue = 0.0;
            double bathymetry = dem_interpolated;

            subcellIndex(0) = xi;
            subcellIndex(1) = yi;

            patch.getAccessor().setValueUNew(subcellIndex, 0, h);
            patch.getAccessor().setValueUNew(subcellIndex, 1, u);
            patch.getAccessor().setValueUNew(subcellIndex, 2, v);
            patch.getAccessor().setParameterWithoutGhostlayer(subcellIndex, 0, mapvalue);
 
            patch.getAccessor().setValueUNew(subcellIndex, 4, h * u);
            patch.getAccessor().setValueUNew(subcellIndex, 5, h * v);
       
            //patch.getAccessor().setParameterWithGhostlayer(subcellIndex, 0, bathymetry);

            // free interpolation
            gsl_interp_free(interp_lower);
            gsl_interp_free(interp_upper);
            gsl_interp_free(interp_temp);

            gsl_interp_accel_free(allocp_lower);
            gsl_interp_accel_free(allocp_upper);
            gsl_interp_accel_free(allocp_temp);
        }
    }

    //Set bathymetry
    update(patch);

//    const tarch::la::Vector<DIMENSIONS, int> subdivisionFactor = patch.getSubdivisionFactor();
//    double min_domainsize = std::min(x_size,y_size);
//    int max_subdivisionFactor = std::max(subdivisionFactor(0),subdivisionFactor(1));
}

tarch::la::Vector<DIMENSIONS,double> peanoclaw::native::scenarios::MekkaFloodSWEScenario::computeDemandedMeshWidth(peanoclaw::Patch& patch, bool isInitializing) {
    double retval = 0.0;

    const tarch::la::Vector<DIMENSIONS, double> patchPosition = patch.getPosition();
    const tarch::la::Vector<DIMENSIONS, double> patchSize = patch.getSize();

    tarch::la::Vector<DIMENSIONS, double> patchCenter;
    patchCenter(0) = patchPosition(0) + patchSize(0)/2.0f;
    patchCenter(1) = patchPosition(1) + patchSize(1)/2.0f;

    const tarch::la::Vector<DIMENSIONS, double> meshWidth = patch.getSubcellSize();
    const tarch::la::Vector<DIMENSIONS, int> subdivisionFactor = patch.getSubdivisionFactor();

    tarch::la::Vector<DIMENSIONS,double> result;
    result = patchSize;
    result[0] = result[0] / (subdivisionFactor(0));
    result[1] = result[1] / (subdivisionFactor(1));
    if (isInitializing) {
        bool refine = check_current_domain(patch);

        std::cout << "using maximal mesh width: " << _maximalMeshWidth[0] << " " << _maximalMeshWidth[1] 
                  << " current patchsize: " << patchSize(0) << " " << patchSize(1)
                  << " current patch position: " << patchPosition(0) << " " << patchPosition(1)
                  << std::endl;

        if (refine) {
            std::cout << "init: we need to refine!" << std::endl;
            result[0] = result[0] / 3;
            result[1] = result[1] / 3;
        } else {
            std::cout << "init: all good" << std::endl;
        }
    } else {
        bool refine = check_current_domain(patch);

        /*if (refine) {
            std::cout << "adaptive: we need to refine!" << std::endl;
            result[0] = result[0] /3;
            result[1] = result[1] /3;
        } else {
            //std::cout << "adaptive: all good" << std::endl;
        }*/
    }

    // control refinement depth

    // minimum depth
    /*result[0] = std::min(_maximalMeshWidth[0], result[0]);
    result[1] = std::min(_maximalMeshWidth[1], result[1]);*/

    // maximum depth
    /*result[0] = std::max(_minialMeshWidth[0], result[0]);
    result[1] = std::max(_minialMeshWidth[1], result[1]);*/
    return result;
}


// TODO: replace loops with standard loops
// TODO: if refinement good enough: fill local patch with values
bool peanoclaw::native::scenarios::MekkaFloodSWEScenario::check_current_domain(
    peanoclaw::Patch& patch
) {
    double min_error = 0.0;
    double max_error = 0.0;
    double maxabs_error = 0.0;

    // --------------------------------------------------- //

    // peano refinement
    size_t refinement[2] = {3,3};

    // TODO: compute pixelspace offset from patch position
    const tarch::la::Vector<DIMENSIONS, double> tree_position = patch.getPosition();
    const tarch::la::Vector<DIMENSIONS, double> patchSize = patch.getSize();

    ssize_t domain_pixelspace_offset[2] = {round(tree_position[0]*scale),round(tree_position[1]*scale)};
  
    const tarch::la::Vector<DIMENSIONS, int> subdivisionFactor = patch.getSubdivisionFactor();
    ssize_t patch_cells[2];
    patch_cells[0] = subdivisionFactor[0];
    patch_cells[1] = subdivisionFactor[1];

    ssize_t patch_vertices[2];
    patch_vertices[0] = patch_cells[0]+1;
    patch_vertices[1] = patch_cells[1]+1;

    double local_pixelspace_size[2];
    local_pixelspace_size[0] = (patchSize[0]*scale) / (patch_vertices[0]-1);
    local_pixelspace_size[1] = (patchSize[1]*scale) / (patch_vertices[1]-1);

    // TODO
    // if (local_pixelspace_size <= 1) abort


    for (ssize_t patch_y = 0; patch_y < patch_cells[1]; ++patch_y) {
        for (ssize_t patch_x = 0; patch_x < patch_cells[0]; ++patch_x) {
            ssize_t position[2] = {patch_x, patch_y};
 
            int pixelspace_x[2];
            int pixelspace_y[2];

            pixelspace_x[0] = domain_pixelspace_offset[0] + round(local_pixelspace_size[0] * position[0]);
            pixelspace_x[1] = domain_pixelspace_offset[0] + round(local_pixelspace_size[0] * (position[0]+1));
            pixelspace_y[0] = domain_pixelspace_offset[1] + round(local_pixelspace_size[1] * position[1]);
            pixelspace_y[1] = domain_pixelspace_offset[1] + round(local_pixelspace_size[1] * (position[1]+1));

            // cell corners
            float dem_height_x0y0 = dem(pixelspace_x[0], pixelspace_y[0]);
            float dem_height_x1y0 = dem(pixelspace_x[1], pixelspace_y[0]);
            float dem_height_x0y1 = dem(pixelspace_x[0], pixelspace_y[1]);
            float dem_height_x1y1 = dem(pixelspace_x[1], pixelspace_y[1]);

#if !defined(NDEBUG)
            std::cout << "position: " << position[0] << " " << position[1] 
                      << " pixelspace x range: " << pixelspace_x[0] << " - " << pixelspace_x[1]-1
                      << " pixelspace y range: " << pixelspace_y[0] << " - " << pixelspace_y[1]-1
                      << " worldspace x range: " << pixelspace_x[0] * dem.scale(0) << " - " << (pixelspace_x[1]-1) * dem.scale(0)
                      << " worldspace y range: " << pixelspace_y[0] * dem.scale(1) << " - " << (pixelspace_y[1]-1) * dem.scale(1)
                      << " heights: " << dem_height_x0y0 << " " << dem_height_x1y0
                      << " " << dem_height_x0y1 << " " << dem_height_x1y1
                      << std::endl;
#endif

            double worldspace_x[2];
            double worldspace_y[2];

            worldspace_x[0] = pixelspace_x[0] * dem.scale(0);
            worldspace_x[1] = pixelspace_x[1] * dem.scale(0);
            worldspace_y[0] = pixelspace_y[0] * dem.scale(1);
            worldspace_y[1] = pixelspace_y[1] * dem.scale(1);

            double lower_heights[2] = {dem_height_x0y0, dem_height_x1y0};
            double upper_heights[2] = {dem_height_x0y1, dem_height_x1y1};

            // setup interpolation
            gsl_interp *interp_lower = gsl_interp_alloc(gsl_interp_linear,2);
            gsl_interp_accel *allocp_lower = gsl_interp_accel_alloc();

            gsl_interp *interp_upper = gsl_interp_alloc(gsl_interp_linear,2);
            gsl_interp_accel *allocp_upper = gsl_interp_accel_alloc();

            gsl_interp_init(interp_lower, worldspace_x, lower_heights, 2);
            gsl_interp_init(interp_upper, worldspace_x, upper_heights, 2);

            gsl_interp *interp_temp = gsl_interp_alloc(gsl_interp_linear,2);
            gsl_interp_accel *allocp_temp = gsl_interp_accel_alloc();


            ssize_t dem_vertices[2];
            dem_vertices[0] = pixelspace_x[1] - pixelspace_x[0];
            dem_vertices[1] = pixelspace_y[1] - pixelspace_y[0];

            // only refine if we actually have some data points left
            if (dem_vertices[0] > refinement[0] && dem_vertices[1] > refinement[1]) {
                for (ssize_t demloop_y=0; demloop_y < dem_vertices[1]; ++demloop_y) {
                    for (ssize_t demloop_x=0; demloop_x < dem_vertices[0]; ++demloop_x) {
                        ssize_t dem_position[2] = {demloop_x, demloop_y};

                        int test_pixelspace_x = dem_position[0] + pixelspace_x[0];
                        int test_pixelspace_y = dem_position[1] + pixelspace_y[0];
                        float dem_height = dem(test_pixelspace_x, test_pixelspace_y);

                        double test_worldspace_x = test_pixelspace_x * dem.scale(0);
                        double test_worldspace_y = test_pixelspace_y * dem.scale(1);

                        double temp[2];
                        temp[0] = gsl_interp_eval(interp_lower, worldspace_x, lower_heights, test_worldspace_x, allocp_lower);
                        temp[1] = gsl_interp_eval(interp_upper, worldspace_x, upper_heights, test_worldspace_x, allocp_upper);

                        gsl_interp_init(interp_temp, worldspace_y, temp, 2);
                        double dem_interpolated = gsl_interp_eval(interp_temp, worldspace_y, temp, test_worldspace_y, allocp_temp);

                        double error = dem_interpolated - dem_height;
                        min_error = std::min(min_error, error);
                        max_error = std::max(max_error, error);
                        maxabs_error = std::max(std::abs(max_error), std::abs(min_error));

                    } // demloop_x
                } // demloop_y
            } else {
#if !defined(NDEBUG)
                std::cout << "reached resolution of dataset" << std::endl;
#endif
            }

            // free interpolation
            gsl_interp_free(interp_lower);
            gsl_interp_free(interp_upper);
            gsl_interp_free(interp_temp);

            gsl_interp_accel_free(allocp_lower);
            gsl_interp_accel_free(allocp_upper);
            gsl_interp_accel_free(allocp_temp);

            // summary
            // go to the next cell
        }
    }

    bool refine = (maxabs_error > 1e1);

#if !defined(NDEBUG)
    if (refine) {
        std::cout << "we need to refine patch: "
                  //<< " min error: " << min_error
                  //<< " max error: " << max_error
                  << " maxabs error: " << maxabs_error
                  << std::endl;
    }
#endif
    return refine;
}
 
void peanoclaw::native::scenarios::MekkaFloodSWEScenario::update(peanoclaw::Patch& patch) {
//    // update bathymetry data
//    //std::cout << "updating bathymetry!" << std::endl;
//    tarch::la::Vector<DIMENSIONS, int> subcellIndex;
//    tarch::la::Vector<DIMENSIONS, double> meshPos;
//    for (int yi = 0; yi < patch.getSubdivisionFactor()(1); yi++) {
//        for (int xi = 0; xi < patch.getSubdivisionFactor()(0); xi++) {
//            subcellIndex(0) = xi;
//            subcellIndex(1) = yi;
//
//            meshPos = patch.getSubcellPosition(subcellIndex);
//            tarch::la::Vector<DIMENSIONS, double> coords = mapMeshToCoordinates(meshPos(0), meshPos(1));
//            double bathymetry = dem((float)coords(0), (float)coords(1));
//            double mapvalue = mapMeshToMap(meshPos);
//
//            patch.setParameterWithGhostlayer(subcellIndex, 0, bathymetry);
//            patch.setParameterWithoutGhostlayer(subcellIndex, 0, mapvalue);
//        }
//    }

#if 0
  tarch::la::Vector<DIMENSIONS, int> subcellIndex;
  tarch::la::Vector<DIMENSIONS, double> meshPos;
  for (int yi = -patch.getGhostlayerWidth(); yi < patch.getSubdivisionFactor()(1) + patch.getGhostlayerWidth(); yi++) {
    for (int xi = -patch.getGhostlayerWidth(); xi < patch.getSubdivisionFactor()(0) + patch.getGhostlayerWidth(); xi++) {
      meshPos = patch.getSubcellPosition(subcellIndex);
      double bathymetry = 0.0;
      if(tarch::la::allGreaterEquals(meshPos, _domainOffset) && !tarch::la::oneGreater(meshPos, _domainOffset+_domainSize)) {
        tarch::la::Vector<DIMENSIONS, double> coords = mapMeshToCoordinates(meshPos(0), meshPos(1));
        bathymetry = dem(coords(0), coords(1));
      }
      subcellIndex(0) = xi;
      subcellIndex(1) = yi;
      patch.getAccessor().setParameterWithGhostlayer(subcellIndex, 0, bathymetry);
    }
  }
#endif

    const tarch::la::Vector<DIMENSIONS, double> tree_position = patch.getPosition();
    const tarch::la::Vector<DIMENSIONS, double> patchSize = patch.getSize();

    ssize_t domain_pixelspace_offset[2] = {round(tree_position[0]*scale),round(tree_position[1]*scale)};
  
    const tarch::la::Vector<DIMENSIONS, int> subdivisionFactor = patch.getSubdivisionFactor();
    ssize_t patch_cells[2];
    patch_cells[0] = subdivisionFactor[0];
    patch_cells[1] = subdivisionFactor[1];

    ssize_t patch_vertices[2];
    patch_vertices[0] = patch_cells[0]+1;
    patch_vertices[1] = patch_cells[1]+1;

    double local_pixelspace_size[2];
    local_pixelspace_size[0] = (patchSize[0]*scale) / (patch_vertices[0]-1);
    local_pixelspace_size[1] = (patchSize[1]*scale) / (patch_vertices[1]-1);


    // initialize new part only
    tarch::la::Vector<DIMENSIONS, int> subcellIndex;
    tarch::la::Vector<DIMENSIONS, double> meshPos;
    std::cout << "update" << std::endl;
    for (int yi = -patch.getGhostlayerWidth(); yi < patch.getSubdivisionFactor()(1) + patch.getGhostlayerWidth(); yi++) {
        for (int xi = -patch.getGhostlayerWidth(); xi < patch.getSubdivisionFactor()(0) + patch.getGhostlayerWidth(); xi++) {
            // evaluate cell centers
            ssize_t position[2] = {xi, yi};
 
            int pixelspace_x[2];
            int pixelspace_y[2];

            pixelspace_x[0] = domain_pixelspace_offset[0] + round(local_pixelspace_size[0] * position[0]);
            pixelspace_x[1] = domain_pixelspace_offset[0] + round(local_pixelspace_size[0] * (position[0]+1));
            pixelspace_y[0] = domain_pixelspace_offset[1] + round(local_pixelspace_size[1] * position[1]);
            pixelspace_y[1] = domain_pixelspace_offset[1] + round(local_pixelspace_size[1] * (position[1]+1));

            // cell corners
            float dem_height_x0y0 = 0.0;
            float dem_height_x1y0 = 0.0;
            float dem_height_x0y1 = 0.0;
            float dem_height_x1y1 = 0.0;
 
            if (pixelspace_x[0] >= 0 && pixelspace_x[1] < dem.dimension(0) 
             && pixelspace_y[0] >= 0 && pixelspace_y[1] <= dem.dimension(1)) 
            {
                dem_height_x0y0 = dem(pixelspace_x[0], pixelspace_y[0]);
                dem_height_x1y0 = dem(pixelspace_x[1], pixelspace_y[0]);
                dem_height_x0y1 = dem(pixelspace_x[0], pixelspace_y[1]);
                dem_height_x1y1 = dem(pixelspace_x[1], pixelspace_y[1]);
            }
     
            double worldspace_x[2];
            double worldspace_y[2];

            worldspace_x[0] = pixelspace_x[0] * dem.scale(0);
            worldspace_x[1] = pixelspace_x[1] * dem.scale(0);
            worldspace_y[0] = pixelspace_y[0] * dem.scale(1);
            worldspace_y[1] = pixelspace_y[1] * dem.scale(1);

            double lower_heights[2] = {dem_height_x0y0, dem_height_x1y0};
            double upper_heights[2] = {dem_height_x0y1, dem_height_x1y1};

            // setup interpolation using cell corners
            gsl_interp *interp_lower = gsl_interp_alloc(gsl_interp_linear,2);
            gsl_interp_accel *allocp_lower = gsl_interp_accel_alloc();

            gsl_interp *interp_upper = gsl_interp_alloc(gsl_interp_linear,2);
            gsl_interp_accel *allocp_upper = gsl_interp_accel_alloc();

            gsl_interp_init(interp_lower, worldspace_x, lower_heights, 2);
            gsl_interp_init(interp_upper, worldspace_x, upper_heights, 2);

            gsl_interp *interp_temp = gsl_interp_alloc(gsl_interp_linear,2);
            gsl_interp_accel *allocp_temp = gsl_interp_accel_alloc();


            ssize_t dem_vertices[2];
            dem_vertices[0] = pixelspace_x[1] - pixelspace_x[0];
            dem_vertices[1] = pixelspace_y[1] - pixelspace_y[0];

            // interpolate on cell center going from lower left corner of cell
            int test_pixelspace_x = pixelspace_x[0];
            int test_pixelspace_y = pixelspace_y[0];

            double test_worldspace_x = (test_pixelspace_x + 0.5 * local_pixelspace_size[0]) * dem.scale(0);
            double test_worldspace_y = (test_pixelspace_y + 0.5 * local_pixelspace_size[1]) * dem.scale(1);

            double temp[2];
            temp[0] = gsl_interp_eval(interp_lower, worldspace_x, lower_heights, test_worldspace_x, allocp_lower);
            temp[1] = gsl_interp_eval(interp_upper, worldspace_x, upper_heights, test_worldspace_x, allocp_upper);

            gsl_interp_init(interp_temp, worldspace_y, temp, 2);
            double dem_interpolated = gsl_interp_eval(interp_temp, worldspace_y, temp, test_worldspace_y, allocp_temp);

            // assign updated bathymetry values
            double bathymetry = dem_interpolated;

            subcellIndex(0) = xi;
            subcellIndex(1) = yi;

            // TODO: value from topological map? into parameter without ghostlayer?
            patch.getAccessor().setParameterWithGhostlayer(subcellIndex, 0, bathymetry);

            // free interpolation
            gsl_interp_free(interp_lower);
            gsl_interp_free(interp_upper);
            gsl_interp_free(interp_temp);

            gsl_interp_accel_free(allocp_lower);
            gsl_interp_accel_free(allocp_upper);
            gsl_interp_accel_free(allocp_temp);
        }
    }
}

tarch::la::Vector<DIMENSIONS,double> peanoclaw::native::scenarios::MekkaFloodSWEScenario::getDomainOffset() const {
  return _domainOffset;
}

tarch::la::Vector<DIMENSIONS,double> peanoclaw::native::scenarios::MekkaFloodSWEScenario::getDomainSize() const {
  return _domainSize;
}

tarch::la::Vector<DIMENSIONS,double> peanoclaw::native::scenarios::MekkaFloodSWEScenario::getInitialMinimalMeshWidth() const {
  return _maximalMeshWidth;
}

tarch::la::Vector<DIMENSIONS,int>    peanoclaw::native::scenarios::MekkaFloodSWEScenario::getSubdivisionFactor() const {
  return _subdivisionFactor;
}

double peanoclaw::native::scenarios::MekkaFloodSWEScenario::getGlobalTimestepSize() const {
  return _globalTimestepSize;
}

double peanoclaw::native::scenarios::MekkaFloodSWEScenario::getEndTime() const {
  return _endTime;
}

double peanoclaw::native::scenarios::MekkaFloodSWEScenario::getInitialTimestepSize() const {
  return 1.0;
}

// box sizes:
// - complete saudi arabia:
//  longitude: 32 - 60E
//  latitude: 13 - 40N
//
// - mekka area


tarch::la::Vector<DIMENSIONS, double> peanoclaw::native::scenarios::MekkaFloodSWEScenario::mapCoordinatesToMesh(double longitude, double latitude) {
    tarch::la::Vector<DIMENSIONS, double> mesh;
 
    /*// put mekkah in the center - adjust bei 0*5 * scale / _domainSize

    double scale = 0.1;

    // peano coordinates
    //float mesh_y = (latitude-13.0f)/27.0f * _domainSize;
    //float mesh_x = (longitude-32.0f)/28.0f * _domainSize;
 
    float mesh_y = (latitude-39.8167f)/1.0f * _domainSize / scale + 0.5;
    float mesh_x = (longitude-21.4167f)/1.0f * _domainSize / scale + 0.5;

    mesh(0) = mesh_x;
    mesh(1) = mesh_y;*/
	

    double lower_left_0 = dem.lower_left(0);
    double lower_left_1 = dem.lower_left(1);
 
    std::cout << "lower left " << lower_left_0 << " " << lower_left_1 << std::endl;

    double ws_0;
    double ws_1;
    dem.pixelspace_to_worldspace(ws_0, ws_1, longitude, latitude);

    // TODO_: the offset change is just due to a problem with plotting, meh ...
    mesh(0) = ws_0 - lower_left_0;
    mesh(1) = ws_1 - lower_left_1;
    return mesh;
}

tarch::la::Vector<DIMENSIONS, double> peanoclaw::native::scenarios::MekkaFloodSWEScenario::mapMeshToCoordinates(double x, double y) {
    tarch::la::Vector<DIMENSIONS, double> coords;

    // put mekkah in the center - adjust bei 0*5 * scale / _domainSize

    /*double scale = 0.5;

    //double latitude = x / _domainSize * 27.0 + 13.0;
    //double longitude = y / _domainSize * 28.0 + 32.0;
 
    double latitude = (y-0.5)* (scale / _domainSize) * 1.0 + 21.4167;
    double longitude = (x-0.5)* (scale / _domainSize) * 1.0 + 39.8167;
    coords(0) = longitude;
    coords(1) = latitude;*/
 
    double lower_left_0 = dem.lower_left(0);
    double lower_left_1 = dem.lower_left(1);
 
    double ps_0 = 0.0;
    double ps_1 = 0.0;
    dem.worldspace_to_pixelspace(ps_0, ps_1, x*scale+lower_left_0, y*scale+lower_left_1);


    coords(0) = ps_0;
    coords(1) = ps_1;
    return coords;
}

double peanoclaw::native::scenarios::MekkaFloodSWEScenario::mapMeshToMap(tarch::la::Vector<DIMENSIONS, double>& coords) {
#if 0
    // relate pixel in png file to bathymetry data
    int width_map = mekka_map.width;
    int height_map = mekka_map.height;
    double width_domain = (dem.upper_right(0) - dem.lower_left(0));
    double height_domain = (dem.upper_right(1) - dem.lower_left(1));

    // relate desired coords to data in map
    double map_pos_x = coords(0)*scale * (width_map / width_domain);
    double map_pos_y = height_map - (coords(1)*scale * (height_map / height_domain)); // picture is upside down otherwise :D

    // now map position to index in png file
    int map_index_x = std::floor(map_pos_x);
    int map_index_y = std::floor(map_pos_y);

    // buffer is organized as grayscale one byte per pixel
    int bufferpos = map_index_y * width_map * 1 + map_index_x * 1;

    // convert gray scale to a double precision value
    double gray_value = mekka_map_data[bufferpos] / 255.0;
#endif
    double gray_value = 0.0;
    return gray_value;
}

#endif
