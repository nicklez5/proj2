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
    state.image_depth = 0;
    state.image_color = new pixel[width*height];
    state.image_depth = new float[width*height];
    for(pixel* xyz = state.image_color; xyz < state.image_color + (width*height); xyz++){
	*xyz = make_pixel(0,0,0);	
    }
    for(int i = 0; i < (width*height);i++){
	state.image_depth[i] = 1;
    }
    
    
   
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
		    for(int xyz = 0 ; xyz < state.num_vertices; xyz += 3){
			clip_triangle(state, (const data_geometry) *state.geometry_array[xyz],(const data_geometry) *state.geometry_array[xyz+1],(const data_geometry) *state.geometry_array[xyz+2],0); 
			//data_geometry *first_triangle = state.geometry_array[xyz];
			//xyz++;
			//data_geometry *second_triangle = state.geometry_array[xyz];
			//xyz++;
			//data_geometry *third_triangle = state.geometry_array[xyz];
			//rasterize_triangle(state,(const data_geometry) *first_triangle,(const data_geometry) *second_triangle ,(const data_geometry) *third_triangle);
				
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
    if(face==1)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    vec4 vertex_a = v0.gl_Position;
    vec4 vertex_b = v1.gl_Position;
    vec4 vertex_c = v2.gl_Position;
    
    data_geometry geo_data[3];
    geo_data[0] = v0;
    geo_data[1] = v1;
    geo_data[2] = v2;
    data_geometry geo_1[3], geo_2[3];
    float a_num1, b_num1, b_num2;
    vec4 p_num1, p_num2;
    if((vertex_a[2] < -vertex_a[3]) && (vertex_b[2] < -vertex_b[3]) && (vertex_c[2] < -vertex_c[3])){
	return;
    }else{
	if((vertex_a[2] < -vertex_a[3]) && (vertex_b[2] >= -vertex_b[3]) && (vertex_c[2] >= -vertex_c[3])){
		b_num1 = (-vertex_b[3] - vertex_b[2]) / (vertex_a[2] + vertex_a[3] - vertex_b[3] - vertex_b[2]);
		b_num2 = (-vertex_a[3] - vertex_a[2]) / (vertex_c[2] + vertex_c[3] - vertex_a[3] - vertex_a[2]);
		p_num1 = b_num1 * vertex_a + (1 - b_num1) * vertex_b;
		p_num2 = b_num2 * vertex_c + (1 - b_num2) * vertex_a;
		geo_1[0].data = new float[state.floats_per_vertex];
		geo_1[1] = v1;
		geo_1[2] = v2;
		for(int i = 0; i < state.floats_per_vertex; i++){
			switch(state.interp_rules[i]){
				case interp_type::flat:
					geo_1[0].data[i] = v0.data[i];
					break;
				case interp_type::smooth:
					geo_1[0].data[i] = b_num2 * v2.data[i] + (1 - b_num2) * v0.data[i];
					break;
				case interp_type::noperspective:
					a_num1 = b_num2 * v2.gl_Position[3] / (b_num2 * v2.gl_Position[3] + (1 - b_num2) * v0.gl_Position[3]);
					geo_1[0].data[i] = a_num1 * v2.data[i] + (1 - a_num1) * v0.data[i];
					break;
				default:
					break;
			}	
		}
		geo_1[0].gl_Position = p_num2;
		geo_data[0] = geo_1[0];
		geo_data[1] = geo_1[1];
		geo_data[2] = geo_1[2];
		clip_triangle(state, geo_data[0],geo_data[1],geo_data[2], face + 1);
		geo_2[0].data = new float[state.floats_per_vertex];
		geo_2[2] = v2;
		
		for(int i = 0; i < state.floats_per_vertex; i++){
			switch(state.interp_rules[i]){
				case interp_type::flat:
					geo_2[0].data[i] = v0.data[i];
					break;
				case interp_type::smooth:
					geo_2[0].data[i] = b_num1 * v0.data[i] + (1-b_num1) * v1.data[i];
					break;
				case interp_type::noperspective:
					a_num1 = b_num1 * v0.gl_Position[3]/(b_num1 * v0.gl_Position[3] + (1-b_num1) * v1.gl_Position[3]);
					geo_2[0].data[i] = a_num1 * v0.data[i] + (1 - a_num1) * v1.data[i];
					break;
				default:
					break;

			}
		}
		geo_2[0].gl_Position = p_num1;
		geo_data[0] = geo_2[0];
		geo_data[1] = geo_1[1];
		geo_data[2] = geo_1[0];

	}
	clip_triangle(state,geo_data[0],geo_data[1],geo_data[2], face + 1);
    }
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    //clip_triangle(state,v0,v1,v2,face+1);
}
float getArea(int x1, int y1, int x2, int y2, int x3, int y3){
	return 0.5 * (
		((x2 * y3 ) - (x3 * y2)) 
		+ ((x3 * y1) - (x1 * y3))
		+ ((x1 * y2) - (y1 * x2))
		);
	
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
    float alpha_p, beta_p, gamma_p;
    int x[3];
    int y[3];
   
    for(int i = 0 ; i < 3 ;i++){
	if(i == 0){
		x[i] = (v0_position[0] * (state.image_width  ) + state.image_width )/2;
		y[i] = (v0_position[1] * (state.image_height ) + state.image_height)/2;
	}else if(i == 1){
		x[i] = (v1_position[0] * (state.image_width ) + state.image_width )/2;
		y[i] = (v1_position[1] * (state.image_height ) + state.image_height )/2;
	}else if(i == 2){
		x[i] = (v2_position[0] * (state.image_width  ) + state.image_width )/2;
		y[i] = (v2_position[1] * (state.image_height  ) + state.image_height )/2;
	}
    }
    int x_min = std::min({x[0],x[1],x[2]});
    int x_max = std::max({x[0],x[1],x[2]});
    int y_min = std::min({y[0],y[1],y[2]});
    int y_max = std::max({y[0],y[1],y[2]});
    x_min = (x_min < 0 ? 0 : x_min);
    x_max = (x_max > state.image_width ? state.image_width  : x_max);
    y_min = (y_min < 0 ? 0 : y_min);
    y_max = (y_max > state.image_height ? state.image_height  : y_max);
   
    
    for(int pixel_x = x_min ; pixel_x < x_max ; pixel_x++){
	for(int pixel_y = y_min; pixel_y < y_max; pixel_y++){
		float alpha_area = getArea(pixel_x, pixel_y, x[1], y[1], x[2], y[2]);
		float beta_area = getArea(x[0],y[0], pixel_x, pixel_y, x[2], y[2]);
		float gamma_area = getArea(x[0],y[0],x[1],y[1],pixel_x,pixel_y);
		//float alpha_area = 0.5*((x[1]*y[2] - x[2]*y[1]) + (x[2]*pixel_y - pixel_x*y[2]) + (pixel_x*y[1] - pixel_y*x[1]));
	       	//float beta_area = 0.5*((pixel_x*y[2] - x[2]*pixel_y) + (x[2]*y[0] - x[0]*y[2]) + (x[0]*pixel_y - pixel_x*y[0]));
		//float gamma_area = 0.5*((x[1]*pixel_y - pixel_x*y[1]) + (pixel_x*y[0] - pixel_y*x[0]) + (x[0]*y[1] - x[1]*y[0]));
		float total_area = alpha_area + beta_area + gamma_area;	
		float _alpha = alpha_area/total_area;
		float _beta = beta_area/total_area;
		float _gamma = gamma_area/total_area;
			
		if((_gamma >= 0 && _gamma <= 1) && (_beta >= 0 && _beta <= 1) && (_alpha >= 0 && _alpha <= 1)){
			data_fragment *temp_fragment = new data_fragment;
			data_output *temp_output = new data_output;
			temp_fragment->data = new float[MAX_FLOATS_PER_VERTEX];
			//temp_fragment->data[0] = pixel_x;
			//temp_fragment->data[1] = pixel_y;
			float depth1 = _alpha * v0_position[2] + _beta * v1_position[2] + _gamma * v2_position[2];
			if(state.image_depth[pixel_y * state.image_width + pixel_x] > depth1){			
				for(int i = 0; i < state.floats_per_vertex ;i++){
					float k_gour;
					switch(state.interp_rules[i]){
						case interp_type::flat:
							temp_fragment->data[i] = v0.data[i];
							//std::cout << "flat rule type" << std::endl;
							break;
						case interp_type::noperspective:
							temp_fragment->data[i] = _alpha * v0.data[i] + _beta * v1.data[i] + _gamma * v2.data[i];							
							//std::cout << "Temp fragment data[" << i << "]: " << temp_fragment->data[i] << std::endl;
							//std::cout << "no perspective rule type" << std::endl;
							break;
							//std::cout << "no perspective rule type" << std::endl;
						case interp_type::smooth:
							k_gour = (_alpha/v0.gl_Position[3] + _beta/v1.gl_Position[3] + _gamma/v2.gl_Position[3]);
							alpha_p = _alpha/(k_gour * v0.gl_Position[3]);
							gamma_p = _gamma/(k_gour* v2.gl_Position[3]);
							beta_p = _beta/(k_gour *v1.gl_Position[3]);
							temp_fragment->data[i] = alpha_p * v0.data[i] + beta_p * v1.data[i] + gamma_p * v2.data[i];
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
				state.image_depth[pixel_y * state.image_width + pixel_x] = depth1;
			}	
		}
			
	
			
			
	}
    }
    
    
    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

