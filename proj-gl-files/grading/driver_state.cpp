#include "driver_state.h"
#include <cstring>
#include <algorithm>
driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color = new pixel[width*height];
    for(pixel* xyz = state.image_color; xyz < state.image_color + (width*height); xyz++){
	*xyz = make_pixel(0,0,0);
    }
    
    
    state.image_depth =0;
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
 
    state.geometry_array = new data_geometry*[state.num_vertices];
    state.vertex_array = new data_vertex*[state.num_vertices];
    for(int i = 0; i < state.num_vertices; i++){
	state.geometry_array[i] = new data_geometry;
	state.vertex_array[i] = new data_vertex;
	state.vertex_array[i]->data = new float[state.floats_per_vertex];
	state.geometry_array[i]->data = new float[state.floats_per_vertex];
    }
    int a = 0;
    for(int i = 0 ; i < state.num_vertices; i++){
	float *data_sub = new float[state.floats_per_vertex];
	for(int sub_a = 0; sub_a < state.floats_per_vertex; sub_a++){
		data_sub[sub_a] = state.vertex_data[a];
		//std::cout << "Print vertex_data: " << state.vertex_data[a] << std::endl;
		a++;
	}
	state.vertex_array[i]->data = data_sub;	
	state.vertex_shader((const data_vertex)*state.vertex_array[i] ,*state.geometry_array[i], state.uniform_data);
        //std::cout << "Geometry vertex_data: " << state.geometry_array[i]->data << std::endl;	
    }
    switch(type){
	    case render_type::triangle:
		    for(int xyz = 0 ; xyz < state.num_vertices; xyz++){
			data_geometry *first_triangle = state.geometry_array[xyz];
			xyz++;
			data_geometry *second_triangle = state.geometry_array[xyz];
			xyz++;
			data_geometry *third_triangle = state.geometry_array[xyz];
			rasterize_triangle(state,(const data_geometry) *first_triangle,(const data_geometry) *second_triangle ,(const data_geometry) *third_triangle);
				
		    }
		    break;
	    case render_type::fan:
		    break;
	    case render_type::strip:
		    break;
	    case render_type::indexed:
		    break;
	    default:
		    break;

    }
    //std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,v0,v1,v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    //std::cout << "Data geometry1 data: " << v0.data << std::endl;
    //std::cout << "Data geometry2 data: " << v1.data << std::endl;
    //std::cout << "Data geometry3 data: " << v2.data << std::endl;
    vec4 v0_position = v0.gl_Position/(v0.gl_Position[3]);
    vec4 v1_position = v1.gl_Position/(v1.gl_Position[3]);
    vec4 v2_position = v2.gl_Position/(v2.gl_Position[3]);

    int x[3];
    int y[3];
    for(int i = 0 ; i < 3 ;i++){
	if(i == 0){
		x[i] = (v0_position[0] * (state.image_width - 1) + state.image_width - 1)/2;
		y[i] = (v0_position[1] * (state.image_height - 1) + state.image_height - 1)/2;
	}else if(i == 1){
		x[i] = (v1_position[0] * (state.image_width - 1) + state.image_width -1)/2;
		y[i] = (v1_position[1] * (state.image_height - 1) + state.image_height - 1)/2;
	}else if(i == 2){
		x[i] = (v2_position[0] * (state.image_width - 1) + state.image_width -1)/2;
		y[i] = (v2_position[1] * (state.image_height - 1) + state.image_height - 1)/2;
	}
    }
    int x_min = std::min({x[0],x[1],x[2]});
    int x_max = std::max({x[0],x[1],x[2]});
    int y_min = std::min({y[0],y[1],y[2]});
    int y_max = std::max({y[0],y[1],y[2]});
    x_min = (x_min < 0 ? 0 : x_min);
    x_max = (x_max > state.image_width ? state.image_width - 1 : x_max);
    y_min = (y_min < 0 ? 0 : y_min);
    y_max = (y_max > state.image_height ? state.image_height - 1 : y_max);
   
    
    for(int pixel_x = x_min ; pixel_x < x_max ; pixel_x++){
	for(int pixel_y = y_min; pixel_y < y_max; pixel_y++){
		float bot_area = 0.5*((x[1]*y[2] - x[2]*y[1]) + (x[2]*y[0] - x[0]*y[2]) + (x[0]*y[1] - x[1]*y[0]));
		float alpha_area = 0.5*((x[1]*y[2] - x[2]*y[1]) + (x[2]*pixel_y - pixel_x*y[2]) + (pixel_x*y[1] - pixel_y*x[1]));
	       	float beta_area = 0.5*((pixel_x*y[2] - x[2]*pixel_y) + (x[2]*y[0] - x[0]*y[2]) + (x[0]*pixel_y - pixel_x*y[0]));
		float gamma_area = 0.5*((x[1]*pixel_y - pixel_x*y[1]) + (pixel_x*y[0] - pixel_y*x[0]) + (x[0]*y[1] - x[1]*y[0]));
		
		float _alpha = alpha_area/bot_area;
		float _beta = beta_area/bot_area;
		float _gamma = gamma_area/bot_area;
			
		if((_gamma >= 0 && _gamma <= 1 ) && (_beta >= 0 && _beta <= 1) && (_alpha >= 0 && _alpha <= 1)){
			data_fragment *temp_fragment = new data_fragment;
			data_output *temp_output = new data_output;
			temp_fragment->data = new float[MAX_FLOATS_PER_VERTEX];
			//temp_fragment->data[0] = pixel_x;
			//temp_fragment->data[1] = pixel_y;
			for(int i = 0; i < state.floats_per_vertex ;i++){
				switch(state.interp_rules[i]){
					case interp_type::flat:
						temp_fragment->data[i] = v0.data[i];
						break;
					case interp_type::noperspective:
						temp_fragment->data[i] = _alpha * v0.data[i] + _beta * v1.data[i] + _gamma * v2.data[i];
						break;
						//std::cout << "no perspective rule type" << std::endl;
						break;
					case interp_type::smooth:
						//std::cout << "Smooth rule type" << std::endl;
						break;
					default:
						//std::cout << "Nothing happened" << std::endl;
						break;
				}
		
			}
			state.fragment_shader(*temp_fragment, *temp_output, state.uniform_data);
			temp_output->output_color = temp_output->output_color *255;
			state.image_color[pixel_y * state.image_width + pixel_x] = make_pixel(temp_output->output_color[0], temp_output->output_color[1], temp_output->output_color[2]);
		
			  
		}
			
			
	}
    }
    
    
    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

