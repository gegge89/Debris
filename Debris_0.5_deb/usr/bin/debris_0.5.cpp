#include <webkit2/webkit2.h>
#include <gtk/gtk.h>
#include <string.h>
#include <gio/gio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream> 
#include <libgen.h>
#include <time.h>
#include <dirent.h>
#include <math.h>
#include <sys/stat.h>
#include <regex>
#include <gst/gst.h>
#include <glib.h>
#include <chrono>
#include <iomanip> 
#include <thread>


#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <netdb.h>
#include <utility>
#include <algorithm>
#include <cctype>
#include <set>
#include <cairo.h>


#include <pwd.h>
#include <cairo/cairo-pdf.h>
#include <gdk-pixbuf/gdk-pixbuf.h>



/*CALCULATOR*/

double expression_result = 0.0;
double first_number = 0.0;
double second_number = 0.0;
int operation = 0; // 0 for none, 1 for addition


GtkWidget *entry_expression;
GtkWidget *entry_result;
GtkWidget *calcmorebutton;
GtkWidget *morebox;
GtkWidget *more_title;

    GtkWidget *second_window;
    GtkWidget *labelparabola_y_formula;
    GtkWidget *labelparabola;
    GtkWidget *label_a_parabola;
    GtkWidget *entrya_parabolay;
    GtkWidget *label_b_parabola;
    GtkWidget *entryb_parabolay;
    GtkWidget *label_c_parabola;
    GtkWidget *entryc_parabolay;
    GtkWidget *label_result_parabolay;
    GtkWidget *label_delta_parabolay;
    GtkWidget *label_vertex_x1_parabolay;
    GtkWidget *label_vertex_x2_parabolay;
    GtkWidget *label_vertex_y_parabolay;
    GtkWidget *label_focus_x_parabolay;
    GtkWidget *label_focus_y_parabolay;
    GtkWidget *sim_axys_parabolay;
    GtkWidget *directrix;
 
    GtkWidget *labelparabola_x_formula;
    GtkWidget *labelparabolax;
    GtkWidget *label_a_parabolax;
    GtkWidget *entrya_parabolax;
    GtkWidget *label_b_parabolax;
    GtkWidget *entryb_parabolax;
    GtkWidget *label_c_parabolax;
    GtkWidget *entryc_parabolax;
    GtkWidget *label_result_parabolax;
    GtkWidget *label_delta_parabolax;
    GtkWidget *label_vertex_x_parabolax;
    GtkWidget *label_vertex_y1_parabolax;
    GtkWidget *label_vertex_y2_parabolax;
    GtkWidget *label_focus_x_parabolax;
    GtkWidget *label_focus_y_parabolax;
    GtkWidget *sim_axys_parabolax;
    GtkWidget *directrix_x;
    
    //GtkWidget *line_box;    
    GtkWidget *label_line_formula;    
    GtkWidget *label_line;  
    GtkWidget *label_m;  
    GtkWidget *label_b;   
    GtkWidget *entry_m;  
    GtkWidget *entry_b; 
    GtkWidget *label_m_line;  
    GtkWidget *label_b_line; 
    GtkWidget *label_result_line; 

    GtkWidget *label_eq_second;  
    GtkWidget *label_eq_second_formula; 
    GtkWidget *label_a_eqsecond;
    GtkWidget *entrya_eqsecond;
    GtkWidget *label_b_eqsecond;
    GtkWidget *entryb_eqsecond;
    GtkWidget *label_c_eqsecond;
    GtkWidget *entryc_eqsecond;
    GtkWidget *label_result_eqsecond; 
    GtkWidget *label_delta_eqsecond; 
    GtkWidget *label_x1_eqsecond; 
    GtkWidget *label_x2_eqsecond; 

    GtkWidget *label_eq_cubic;  
    GtkWidget *label_eq_cubic_formula; 
    GtkWidget *entrya_eqcubic;
    GtkWidget *entryb_eqcubic;
    GtkWidget *entryc_eqcubic;
    GtkWidget *entryd_eqcubic;
    GtkWidget *label_result_eqcubic; 
    GtkWidget *label_x_eqcubic;  

    GtkWidget *label_eq_archim_formula; 
    GtkWidget *entryt_eqarchim;
    GtkWidget *entryv_eqarchim;
    GtkWidget *entryw_eqarchim;
    
    GtkWidget *label_expo_formula;
    GtkWidget *entryx1_expo;
    GtkWidget *entryx2_expo;
    GtkWidget *entryx3_expo;
    GtkWidget *entryx4_expo;
    GtkWidget *entryx5_expo;
    GtkWidget *entryx6_expo;  
    GtkWidget *entryb_expo;    

    GtkWidget *label_a_proportions;
    GtkWidget *label_b_proportions;
    GtkWidget *label_c_proportions;  
    GtkWidget *label_d_proportions;  
    GtkWidget *entrya_proportions;
    GtkWidget *entryb_proportions;
    GtkWidget *entryc_proportions;  
    GtkWidget *entryd_proportions;  
    GtkWidget *label_x_proportions;
    GtkWidget *label_result_proportions;
    GtkWidget *entries_proportions;

    GtkWidget *label_a_ineq2;
    GtkWidget *label_b_ineq2;
    GtkWidget *label_c_ineq2;  
    GtkWidget *label_x_ineq2;
    GtkWidget *label_y_ineq2;
    GtkWidget *label_x2_ineq2;
    GtkWidget *label_y2_ineq2;
    GtkWidget *label_result_ineq2;

 GtkEntry *number;
    
static gboolean draw_cos_wave = FALSE;
static gboolean draw_cos_wave2 = FALSE;
static gdouble calculated_cosine = 0.0;

static gboolean draw_sin_wave = FALSE;
static gboolean draw_sin_wave2 = FALSE;
static gboolean draw_sin_wave3 = FALSE;
static gboolean draw_sin_wave4 = FALSE;
static gdouble calculated_sin = 0.0;

static gboolean draw_atan = FALSE;
static gdouble calculated_atan = 0.0;
static gboolean draw_tan = FALSE;
static gdouble calculated_tan = 0.0;
static gdouble origin = 0.0;

static gboolean draw_parabolay = FALSE;
static gboolean draw_parabolax = FALSE;
static gdouble py_delta = 0.0;
static gdouble py_xvertex = 0.0;
static gdouble py_x2vertex = 0.0;
static gdouble py_yvertex = 0.0;
static gdouble py_y2vertex = 0.0;
static gdouble py_xfocus = 0.0;
static gdouble py_yfocus = 0.0;
static gdouble py_simaxys = 0.0;
static gdouble py_dir = 0.0;

static gboolean draw_line = FALSE;
static gdouble res_m_line = 0.0;
static gdouble res_b_line = 0.0;

static gboolean draw_eq2 = FALSE;
static gdouble a_eq2 = 0.0; //a
static gdouble b_eq2 = 0.0; //b
static gdouble c_eq2 = 0.0; //c
static gdouble x1_eq2 = 0.0; //x1
static gdouble x2_eq2 = 0.0; //x2
static gdouble delta_eq2 = 0.0; //delta

static gboolean draw_eq3 = FALSE;
static gdouble a_eq3 = 0.0; //a
static gdouble b_eq3 = 0.0; //b
static gdouble c_eq3 = 0.0; //c
static gdouble d_eq3 = 0.0; //c
static gdouble x_eq3 = 0.0; //x1

static gboolean draw_archim = FALSE;
static gdouble t_eqarchim = 0.0; //a
static gdouble v_eqarchim = 0.0; //b
static gdouble w_eqarchim = 0.0; //c

static gboolean draw_expo = FALSE;
static gdouble x1_expo = 0.0; 
static gdouble x2_expo = 0.0; 
static gdouble x3_expo = 0.0; 
static gdouble x4_expo = 0.0; 
static gdouble x5_expo = 0.0; 
static gdouble x6_expo = 0.0; 
static gdouble b_expo = 0.0; 

static gboolean draw_proportions = FALSE;
static gdouble a_proportions = 0.0; 
static gdouble b_proportions = 0.0; 
static gdouble c_proportions = 0.0; 
static gdouble d_proportions = 0.0; 
static gdouble x_proportions = 0.0; 

static gboolean draw_ineq2 = FALSE;
static gdouble a_ineq2 = 0.0; 
static gdouble b_ineq2 = 0.0; 
static gdouble c_ineq2 = 0.0; 
static gdouble x1_ineq2 = 0.0; 
static gdouble x2_ineq2 = 0.0; 
static gdouble y1_ineq2 = 0.0; 
static gdouble y2_ineq2 = 0.0; 


void drawResultCircle(cairo_t *cr, bool draw_circle, double expression_result, double width, double height, double scale_factor, double step) {
    if (!draw_circle) {
        return; // Skip drawing the circle if the flag is false
    }

    // Set the color to red
    cairo_set_source_rgb(cr, 1, 0, 0);

    // Calculate the position
    double x_position = width / 2 + (expression_result / scale_factor) * step;
    double y_position = height / 2 - (expression_result / scale_factor) * step;

    // Ensure the circle is within bounds before drawing
    if (x_position >= 0 && x_position <= width && y_position >= 0 && y_position <= height) {
        cairo_arc(cr, x_position, y_position, 5, 0, 2 * G_PI); // Draw the circle
        cairo_fill(cr); // Fill the circle
    }
}

static void draw_callback(GtkWidget *widget, cairo_t *cr, gpointer data) {
    // Get widget dimensions
    int width = gtk_widget_get_allocated_width(widget);
    int height = gtk_widget_get_allocated_height(widget);

    // Set background color
    cairo_set_source_rgb(cr, 38.0 / 255, 38.0 / 255, 38.0 / 255); // #262626
   // cairo_set_source_rgb(cr, 38.0 / 189, 38.0 / 189, 38.0 / 189); 
    cairo_paint(cr);

// Define the offsets
int vertical_offset = -10;  // Adjust this value to move the grid up or down
int horizontal_offset = -3; // Adjust this value to move the grid left or right

// Draw grid lines
cairo_set_source_rgb(cr, 0, 0, 0); // Black
cairo_set_line_width(cr, 1);

for (int i = 0; i < width; i += 20) {
    cairo_move_to(cr, i + horizontal_offset, vertical_offset);
    cairo_line_to(cr, i + horizontal_offset, height + vertical_offset);
}

for (int i = 0; i < height; i += 20) {
    cairo_move_to(cr, horizontal_offset, i + vertical_offset);
    cairo_line_to(cr, width + horizontal_offset, i + vertical_offset);
}
cairo_stroke(cr);


    // Draw axes
    cairo_set_source_rgb(cr, 0.0, 170.0 / 255, 255.0 / 255); // #00aaff
    cairo_set_line_width(cr, 2);

    // X-axis
    cairo_move_to(cr, 0, height / 2);
    cairo_line_to(cr, width, height / 2);
    cairo_stroke(cr);

    // Y-axis
    cairo_move_to(cr, width / 2, 0);
    cairo_line_to(cr, width / 2, height);
    cairo_stroke(cr);

    // Draw axis labels
    cairo_set_source_rgb(cr, 1, 1, 1); // White
    cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 15);

    cairo_move_to(cr, width - 20, height / 2 - 10);
    cairo_show_text(cr, "X");

    cairo_move_to(cr, width / 2 + 10, 20);
    cairo_show_text(cr, "Y");

    // Determine scale factor
    int scale_factor = 1;
    if (expression_result > 9 && expression_result <= 90) {
        scale_factor = 10;
    } else if (expression_result > 90 && expression_result <= 200) {
        scale_factor = 20;
    } else if (expression_result > 200 && expression_result <= 499) {
        scale_factor = 50;
    } else if (expression_result > 499 && expression_result <= 999) {
        scale_factor = 100;
    } else if (expression_result > 999 && expression_result <= 9999) {
        scale_factor = 1000;
    } else if (expression_result > 9999 && expression_result <= 99999) {
        scale_factor = 10000;
    } else if (expression_result > 99999) {
        scale_factor = 100000;
    }

    // Declare step variable
    int step = 20;

    // Draw numbers along x and y axes
    if (!(draw_cos_wave || draw_cos_wave2 || draw_sin_wave || draw_sin_wave2 || draw_sin_wave3 || draw_sin_wave4)) {
        cairo_set_font_size(cr, 12);  
        for (int i = -width / 2 / step; i <= width / 2 / step; ++i) {
            if (i != 0) {
                char label[10];
                snprintf(label, sizeof(label), "%d", i * scale_factor);
                cairo_move_to(cr, width / 2 + i * step - 5, height / 2 + 15);
                cairo_show_text(cr, label);
            }
        }

        for (int i = -height / 2 / step; i <= height / 2 / step; ++i) {
            if (i != 0) {
                char label[10];
                snprintf(label, sizeof(label), "%d", -i * scale_factor);
                cairo_move_to(cr, width / 2 + 5, height / 2 + i * step + 5);
                cairo_show_text(cr, label);
            }
        }

        // Draw a circle representing the expression result
        drawResultCircle(cr, true, expression_result, width, height, scale_factor, step);
    }

    if (draw_archim) {
        cairo_save(cr); // Save the current state

        // Center of the Cartesian plane
        double centerX = width / 2;
        double centerY = height / 2;

        // Time, velocity, and angular velocity
        double t = t_eqarchim;
        double v = v_eqarchim;
        double w = w_eqarchim;

        // Calculate the constant k
        double k = v / w;

        // Set the color for the spiral
        cairo_set_source_rgb(cr, 1.0, 0.0, 1.0); //purple but  Orange

        // Draw the spiral
        cairo_move_to(cr, centerX, centerY);
        for (double theta = 0; theta < t * w; theta += 0.01) {
            double r = k * theta;
            double x = centerX + r * cos(theta);
            double y = centerY + r * sin(theta);
            cairo_line_to(cr, x, y);
        }

        cairo_stroke(cr);
        cairo_restore(cr); // Restore the saved state
    }

if (draw_line) {
    cairo_save(cr); // Save the current state

    //g_print("function is called LINE\n");

    cairo_set_source_rgb(cr, 1.0, 0.0, 1.0); // purple
    cairo_set_line_width(cr, 2.0); // Set line width to a reasonable value, e.g., 2.0

    // Calculate the coordinates for the line based on the values of res_m_line and res_b_line
   double x1, y1, x2, y2;
    if (res_m_line > 0) {
        if (res_b_line > 0) {
            x1 = width / 2 + (-(res_b_line / res_m_line) ) * 20.0;
        } else if (res_b_line < 0) {
            x1 = width / 2 + (-(res_b_line + res_m_line) ) * 20.0;
        } else {
            x1 = width / 2; // Default value if res_b_line is 0
        }
        y1 = height / 2;
        x2 = width / 2 + 0 * 20.0; // Point two x-coordinate, scaled (0 remains 0)
        y2 = height / 2 - res_b_line * 20.0; // Point two y-coordinate, scaled
    } else if (res_m_line < 0) {
        if (res_b_line > 0) {
            x1 = width / 2 + (-(res_b_line / res_m_line)) * 20.0;
        } else if (res_b_line < 0) {
            x1 = width / 2 + (-(res_b_line / res_m_line)) * 20.0;
        } else {
            x1 = width / 2; // Default value if res_b_line is 0
        }
        y1 = height / 2;
        x2 = width / 2 + 0 * 20.0; // Point two x-coordinate, scaled (0 remains 0)
        y2 = height / 2 - res_b_line * 20.0; // Point two y-coordinate, scaled
    } else { // res_m_line == 0
        x1 = 0;
        y1 = height / 2 - res_b_line * 20.0; // y-coordinate is res_b_line
        x2 = width;
        y2 = height / 2 - res_b_line * 20.0; // y-coordinate is res_b_line
    }

//START LINE
  if (res_m_line != 0) {
        // Calculate the slope (m) of the line
        double dx = x2 - x1;
        double dy = y2 - y1;
        double m = dy / dx;

        double x_start = 0;
        double y_start = y1 + m * (x_start - x1);
        double x_end = width;
        double y_end = y1 + m * (x_end - x1);

        // Move to the starting point (x_start, y_start)
        cairo_move_to(cr, x_start, y_start);

        // Draw the line to the ending point (x_end, y_end)
        cairo_line_to(cr, x_end, y_end);
    } else {
        // Draw a straight line parallel to the x-axis
        cairo_move_to(cr, x1, y1);
        cairo_line_to(cr, x2, y2);
    }

    // Stroke the line
    cairo_stroke(cr);

    cairo_restore(cr); // Restore the saved state
}
//END LINE

//START PARABOLA Y
if (draw_parabolay) {
    cairo_save(cr); // Save the current state

    // Circle color and radius
    double circle_radius = 5.0; // Radius of the circle

    // Circles are already correctly placed
    cairo_set_source_rgba(cr, 0, 0, 1, 1.0); // Blue for the focus point
    cairo_arc(cr, width / 2 + py_xfocus * 20.0, height / 2 - py_yvertex * 20.0, circle_radius, 0, 2 * M_PI);
    cairo_fill(cr);

    cairo_set_source_rgba(cr, 0, 1, 0, 1.0); // Green for the first root
    cairo_arc(cr, width / 2 + py_xvertex * 20.0, height / 2, circle_radius, 0, 2 * M_PI); // y = 0
    cairo_fill(cr);

    cairo_set_source_rgba(cr, 1, 0, 0, 1.0); // Red for the second root
    cairo_arc(cr, width / 2 + py_x2vertex * 20.0, height / 2, circle_radius, 0, 2 * M_PI); // y = 0
    cairo_fill(cr);

    // Use coordinates for key points to dynamically calculate the parabola
    double x1 = py_xvertex;       // First root
    double x2 = py_x2vertex;      // Second root
    double x_vertex = py_xfocus;  // x-coordinate of the vertex
    double y_vertex = py_yvertex; // y-coordinate of the vertex

    // Calculate the parabola coefficient 'a' using vertex and roots
    double a = y_vertex / ((x_vertex - x1) * (x_vertex - x2));

    // Draw the parabola curve
    cairo_set_source_rgba(cr, 0.5, 0, 0.5, 1.0); // Purple color for the parabola
    cairo_set_line_width(cr, 3);

    // Start drawing at x = -12
    cairo_move_to(cr, width / 2 + (-12) * 20.0, height / 2 - (a * (-12 - x1) * (-12 - x2)) * 20.0);

    for (double x = -12; x <= 12; x += 0.01) {
        double y = a * (x - x1) * (x - x2); // Parabola equation
        cairo_line_to(cr, width / 2 + x * 20.0, height / 2 - y * 20.0);
    }

    cairo_stroke(cr); // Render the parabola curve

    cairo_restore(cr); // Restore the saved state
}

//END PARABOLA Y

//START PARABOLA X
if (draw_parabolax) {
    cairo_save(cr); // Save the current state

    // Circle color and radius
    double circle_radius = 5.0; // Radius for the circles

    // Draw blue circle for the vertex/focus at (py_xvertex, py_xfocus)
    cairo_set_source_rgba(cr, 0, 0, 1, 1.0); // Blue color
    cairo_arc(cr, width / 2 + py_xvertex * 20.0,
                 height / 2 - py_xfocus * 20.0,
                 circle_radius, 0, 2 * M_PI);
    cairo_fill(cr);

    // Draw green circle for the first root (at x = 0, y = py_yvertex)
    cairo_set_source_rgba(cr, 0, 1, 0, 1.0); // Green color
    cairo_arc(cr, width / 2,
                 height / 2 - py_yvertex * 20.0,
                 circle_radius, 0, 2 * M_PI);
    cairo_fill(cr);

    // Draw red circle for the second root (at x = 0, y = py_y2vertex)
    cairo_set_source_rgba(cr, 1, 0, 0, 1.0); // Red color
    cairo_arc(cr, width / 2,
                 height / 2 - py_y2vertex * 20.0,
                 circle_radius, 0, 2 * M_PI);
    cairo_fill(cr);

    // Set up key points for the parabola.
    // For a parabola where x = a*(y - y1)*(y - y2), the roots are at y1 and y2,
    // and we require that when y equals the vertex’s y,
    // x becomes py_xvertex. Here the vertex is defined at (py_xvertex, py_xfocus).
    double y1 = py_yvertex;       // First root (y-axis coordinate)
    double y2 = py_y2vertex;      // Second root (y-axis coordinate)
    double y_vertex = py_xfocus;  // Use py_xfocus so that the vertex is at the blue circle
    double x_vertex = py_xvertex; // x-coordinate of the vertex

    // Calculate coefficient 'a' so that x = a*(y - y1)*(y - y2)
    // gives x = py_xvertex when y = py_xfocus.
    double a = x_vertex / ((y_vertex - y1) * (y_vertex - y2));

    // Set up drawing for the parabola curve
    cairo_set_source_rgba(cr, 0.5, 0, 0.5, 1.0); // Purple for the parabola
    cairo_set_line_width(cr, 3);

    // Start drawing the curve at y = -12 (just as in your other function)
    cairo_move_to(cr, width / 2 + (-12) * 20.0,
                      height / 2 - (a * (-12 - y1) * (-12 - y2)) * 20.0);

    // Draw the parabola for y in [-12, 12]
    for (double y = -12; y <= 12; y += 0.01) {
        double x = a * (y - y1) * (y - y2);  // Parabola equation: x = a*(y - y1)*(y - y2)
        cairo_line_to(cr, width / 2 + x * 20.0,
                          height / 2 - y * 20.0);
    }

    cairo_stroke(cr); // Render the parabola

    cairo_restore(cr); // Restore the previous state
}


//END PARABOLA X

////SQUARE EQUATION
if (draw_eq2) {
    cairo_save(cr); // Save the current state
    
    // Calculate vertex coordinates without the translation
    double x_vertex = -b_eq2 / (2 * a_eq2);
    double y_vertex = a_eq2 * pow(x_vertex, 2) + b_eq2 * x_vertex + c_eq2;
    
    // Calculate the x-intercepts (roots of the equation)
    double x1 = x1_eq2;
    double x2 = x2_eq2;
    
    //g_print("x1: %f\n", x1);
    //g_print("x2: %f\n", x2);
    
    cairo_set_source_rgb(cr, 1.0, 0.0, 1.0); // Orange

    // Draw circles at the x-intercepts
    cairo_arc(cr, width / 2 + x1 * 20.0, height / 2 - 0.0 * 20.0, 5, 0, 2 * M_PI); 
    cairo_fill(cr);

    cairo_arc(cr, width / 2 + x2 * 20.0, height / 2 - 0.0 * 20.0, 5, 0, 2 * M_PI); 
    cairo_fill(cr);

    cairo_set_source_rgba(cr, 1, 0, 1, 1.0); // Purple for the parabola curve with full opacity
    cairo_set_line_width(cr, 3);

    // Draw the parabola
    cairo_move_to(cr, width / 2 + (-12) * 20.0, height / 2 - (a_eq2 * pow(-12, 2) + b_eq2 * -12 + c_eq2) * 20.0); 

    for (double x = -12; x <= 12; x += 0.01) {
        double y = a_eq2 * pow(x, 2) + b_eq2 * x + c_eq2;
        cairo_line_to(cr, width / 2 + x * 20.0, height / 2 - y * 20.0);
    }

    cairo_stroke(cr);
    cairo_restore(cr); // Restore the saved state          
}

/////SQAURE EQUATION END

//CUBIC EQUATION
    if (draw_eq3) {
        cairo_save(cr); // Save the current state
        
        // Function to draw cubic equation
        // Coefficients a_eq3, b_eq3, c_eq3, d_eq3 defined elsewhere in the code
        double a = a_eq3;
        double b = b_eq3;
        double c = c_eq3;
        double d = d_eq3;

        // Define scale and offset
        double scale_x = step;
        double scale_y = step;
        double offset_x = width / 2;
        double offset_y = height / 2;

        cairo_set_source_rgb(cr, 1.0, 0.0, 1.0); 
        cairo_set_line_width(cr, 2);
        
        for (double x = -offset_x / scale_x; x <= offset_x / scale_x; x += 0.01) {
            double y = a * x * x * x + b * x * x + c * x + d;
            double draw_x = offset_x + x * scale_x;
            double draw_y = offset_y - y * scale_y;
            if (x == -offset_x / scale_x) {
                cairo_move_to(cr, draw_x, draw_y);
            } else {
                cairo_line_to(cr, draw_x, draw_y);
            }
        }
        cairo_stroke(cr);

        cairo_restore(cr); // Restore the saved state
    }


//END CUBIC EQUATION
//START EXPONENTIAL FUNCTION
if (draw_expo) {
    cairo_save(cr); // Save the current state

    // Set the color for the exponential curve.
    // (Currently magenta; adjust if you prefer green.)
    cairo_set_source_rgb(cr, 1.0, 0.0, 1.0);
    cairo_set_line_width(cr, 2.0);

    // Array of 6 x-values where circles will be drawn.
    gdouble x_values[6] = {x1_expo, x2_expo, x3_expo, x4_expo, x5_expo, x6_expo};

    // Define conversion factors (pixels per unit) and the origin of the Cartesian plane.
    const double scale = 20.0;
    const double origin_x = width / 2.0;
    const double origin_y = height / 2.0;

    // Draw the exponential curve f(x) = b^x.
    // Use a step size for a smoother curve.
    double step_size = 0.1;
    bool first_point = true; // flag for move_to instead of line_to

    for (double x = -10.0; x < 10.0; x += step_size) {
        double y = pow(b_expo, x);
        double x_pixels = origin_x + x * scale;
        double y_pixels = origin_y - y * scale; // Invert the y-axis

        if (first_point) {
            cairo_move_to(cr, x_pixels, y_pixels);
            first_point = false;
        } else {
            cairo_line_to(cr, x_pixels, y_pixels);
        }
    }
    cairo_stroke(cr);

    // Draw circles at each defined x-value.
    // These circles mark key points on the curve.
    cairo_set_source_rgba(cr, 0, 1, 0, 1.0); // Green color for circles
    cairo_set_line_width(cr, 3);

    for (int i = 0; i < 6; ++i) {
        double x = x_values[i];
        double y = pow(b_expo, x);
        double x_pixels = origin_x + x * scale;
        double y_pixels = origin_y - y * scale;

        // Draw a circle with a radius of 3 pixels.
        cairo_arc(cr, x_pixels, y_pixels, 3, 0, 2 * G_PI);
        cairo_fill(cr);
    }

    cairo_restore(cr); // Restore the previous state
}


//END EXPONENTIAL FUNCTION
 
//START PROPORTIONS 
    if (draw_proportions) {
cairo_save(cr);

// Function call where the circle should NOT be drawn
drawResultCircle(cr, false, expression_result, width, height, scale_factor, step);

        // Circle one
        cairo_set_source_rgba(cr, 0, 1, 0, 1.0); // Green for the position with full opacity
        double x_coord1 = width / 2 + (a_proportions / scale_factor) * step;
        double y_coord1 = height / 2 - (b_proportions / scale_factor) * step;
        double radius = 5.0; // Radius of the circle
        if (x_coord1 >= 0 && x_coord1 <= width && y_coord1 >= 0 && y_coord1 <= height) {
            cairo_arc(cr, x_coord1, y_coord1, radius, 0, 2 * M_PI); // Draw the circle
            cairo_fill(cr); // Fill the circle

            // Add label with values
            
            gchar *label1 = g_strdup_printf("(%.2f, %.2f)", a_proportions, b_proportions);
            cairo_set_source_rgb(cr, 0, 1, 0); // green text
            cairo_move_to(cr, width / 2 + 50, height / 2 + 80);
            cairo_show_text(cr, label1);
            g_free(label1);
        }

        // Circle two
        cairo_set_source_rgba(cr, 1, 0, 1, 1.0); // Red for the position with full opacity
        double x_coord2 = width / 2 + (c_proportions / scale_factor) * step;
        double y_coord2 = height / 2 - (d_proportions / scale_factor) * step;
        if (x_coord2 >= 0 && x_coord2 <= width && y_coord2 >= 0 && y_coord2 <= height) {
            cairo_arc(cr, x_coord2, y_coord2, radius, 0, 2 * M_PI); // Draw the circle
            cairo_fill(cr); // Fill the circle

            // Add label with values
            
            gchar *label2 = g_strdup_printf("(%.2f, %.2f)", c_proportions, d_proportions);
            cairo_set_source_rgb(cr, 1, 0, 1); 
            cairo_move_to(cr, width / 2 + 50, height / 2 + 100); 
            cairo_show_text(cr, label2);
            g_free(label2);
        }

            gchar *label3 = g_strdup_printf("Hightest value");
            cairo_set_source_rgb(cr, 1, 0, 0); 
            cairo_move_to(cr, width / 2 + 50, height / 2 + 120); 
            cairo_show_text(cr, label3);
            g_free(label3);

        cairo_stroke(cr);
        cairo_save(cr);
        cairo_restore(cr); // Restore the saved state
    }

//END PROPORTIONS

//START INEQ2
    if (draw_ineq2) {
        cairo_save(cr); 

        // Calculate the coordinates for drawing
        double x_coord1 = width / 2 + (x1_ineq2 / scale_factor) * step;
        double y_coord1 = height / 2 - (y1_ineq2 / scale_factor) * step;
        double x_coord2 = width / 2 + (x2_ineq2 / scale_factor) * step;
        double y_coord2 = height / 2 - (y2_ineq2 / scale_factor) * step;

        // Calculate the slope of the line
        double slope = (y_coord2 - y_coord1) / (x_coord2 - x_coord1);

        // Calculate the points where the line intersects the canvas edges
        double x_left = 0;
        double y_left = y_coord1 - (slope * (x_coord1 - x_left));
        
        double x_right = width;
        double y_right = y_coord1 + (slope * (x_right - x_coord1));
        
        double y_top = 0;
        double x_top = x_coord1 - (y_coord1 - y_top) / slope;
        
        double y_bottom = height;
        double x_bottom = x_coord1 + (y_bottom - y_coord1) / slope;

        // Draw the shaded area above the line
        cairo_set_source_rgba(cr, 0, 0, 1, 0.3); // Blue with low opacity
        cairo_move_to(cr, x_left, y_left);
        cairo_line_to(cr, x_right, y_right);
        cairo_line_to(cr, x_right, 0);
        cairo_line_to(cr, 0, 0);
        cairo_close_path(cr);
        cairo_fill(cr);

        // Draw the shaded area below the line
        cairo_set_source_rgba(cr, 1, 0, 0, 0.3); // Red with low opacity
        cairo_move_to(cr, x_left, y_left);
        cairo_line_to(cr, x_right, y_right);
        cairo_line_to(cr, x_right, height);
        cairo_line_to(cr, 0, height);
        cairo_close_path(cr);
        cairo_fill(cr);

        // Draw the circles
        cairo_set_source_rgba(cr, 0, 1, 0, 1.0); // Green for the first circle
        double radius = 5.0; // Radius of the circle
        if (x_coord1 >= 0 && x_coord1 <= width && y_coord1 >= 0 && y_coord1 <= height) {
            cairo_arc(cr, x_coord1, y_coord1, radius, 0, 2 * M_PI); // Draw the circle
            cairo_fill(cr); // Fill the circle

            // Add label with values
            gchar *label1 = g_strdup_printf("(%.2f, %.2f)", x1_ineq2, y1_ineq2);
            cairo_set_source_rgb(cr, 1, 1, 1); // White for text
            cairo_move_to(cr, x_coord1 + radius + 2, y_coord1 - radius - 2);
            cairo_show_text(cr, label1);
            g_free(label1);
        }

        cairo_set_source_rgba(cr, 1, 0, 0, 1.0); // Red for the second circle
        if (x_coord2 >= 0 && x_coord2 <= width && y_coord2 >= 0 && y_coord2 <= height) {
            cairo_arc(cr, x_coord2, y_coord2, radius, 0, 2 * M_PI); // Draw the circle
            cairo_fill(cr); // Fill the circle

            // Add label with values
            gchar *label2 = g_strdup_printf("(%.2f, %.2f)", x2_ineq2, y2_ineq2);
            cairo_set_source_rgb(cr, 1, 1, 1); // White for text
            cairo_move_to(cr, x_coord2 + radius + 2, y_coord2 - radius - 2);
            cairo_show_text(cr, label2);
            g_free(label2);
        }

        // Draw the infinite line
        cairo_set_source_rgb(cr, 1, 1, 1); // white for the line

        // Determine which points to use for the infinite line
        if (y_left >= 0 && y_left <= height) {
            cairo_move_to(cr, x_left, y_left);
        } else {
            cairo_move_to(cr, x_top, y_top);
        }

        if (y_right >= 0 && y_right <= height) {
            cairo_line_to(cr, x_right, y_right);
        } else {
            cairo_line_to(cr, x_bottom, y_bottom);
        }

        cairo_stroke(cr); // Draw the line

        cairo_restore(cr); // Restore the saved state
    }


//END INEQ2

//START COS&SIN
    // Additional drawing logic when specific waves are active
    if (draw_cos_wave || draw_cos_wave2 || draw_sin_wave || draw_sin_wave2 || draw_sin_wave3 || draw_sin_wave4) {
        cairo_save(cr); // Save the current state

        cairo_translate(cr, width / 2, height / 2);
        cairo_scale(cr, height / 2.0, height / 2.0); // Scale uniformly based on the smaller dimension

        cairo_new_path(cr);
        // Draw a circle with radius 1
        cairo_set_source_rgb(cr, 0.0, 0.6667, 1.0); // Blue
        cairo_set_line_width(cr, 0.02); // Adjust line width for better visibility
        // Draw the circle
        cairo_arc(cr, 0, 0, 1, 0, 2 * M_PI); // Complete the circle
        cairo_stroke(cr); // Stroke the circle

        // Draw numbers along x and y axes
        cairo_set_source_rgb(cr, 1.0, 1.0, 1.0); // White
        cairo_set_font_size(cr, 0.05); // Smaller font size

        for (double i = -1.0; i <= 1.0; i += 0.5) {
            char label[10];
            snprintf(label, sizeof(label), "%.1f", i);
            cairo_move_to(cr, i - 0.05, 0.1);
            cairo_show_text(cr, label);
            cairo_move_to(cr, 0.05, -i + 0.05);
            cairo_show_text(cr, label);
        }

        // Calculate cosine and sine for the current expression result
        double cos_result = cos(expression_result * M_PI / 180.0);
        double sin_result = sin(expression_result * M_PI / 180.0);

        double cos_x_position = expression_result; // Cosine for x-axis
        double sin_y_position = expression_result * -1; // Sine for y-axis

        // Draw the cosine wave position
if (draw_cos_wave) {
    if (cos_x_position >= -1 && cos_x_position <= 1) {
        // Draw the cosine result dot
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the position
        cairo_arc(cr, cos_x_position, 0, 0.02, 0, 2 * M_PI); // Draw a small circle at the position
        cairo_fill(cr);

        // Draw a vertical line from the cosine result dot to the origin (0, 0)
        cairo_set_source_rgb(cr, 1, 0, 1); // Purple for the vertical line
        cairo_set_line_width(cr, 0.01); // Set line width to a thinner value
        cairo_move_to(cr, cos_x_position, 0); // Start from the dot position
        cairo_line_to(cr, 0, 0); // Draw the line to the origin (0, 0)
        cairo_stroke(cr);

        // Calculate the angle corresponding to the cosine value
        double angle = acos(cos_x_position); // Angle in radians (0 to π)
        angle = 2 * M_PI - angle; // Adjust angle for correct placement (360° - angle)
        double circle_x = cos(angle); // X-coordinate on the unit circle
        double circle_y = sin(angle); // Y-coordinate on the unit circle

        // Draw a line from the point on the circle to the center
        cairo_set_source_rgb(cr, 0, 0, 1); // Blue for the circle-to-center line
        cairo_set_line_width(cr, 0.01); // Maintain the same line width
        cairo_move_to(cr, circle_x, circle_y); // Start from the calculated point on the circle
        cairo_line_to(cr, 0, 0); // Draw to the center
        cairo_stroke(cr);
        
        // Draw the green line from the cosine dot to the circle point
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the line
        cairo_set_line_width(cr, 0.01); // Maintain thinner line width
        cairo_move_to(cr, cos_x_position, 0); // Start from the dot position
        cairo_line_to(cr, circle_x, circle_y); // Stop at the circle point
        cairo_stroke(cr);       
        
//labels        
gchar *cos_label = g_strdup_printf("Cosine");
cairo_set_source_rgb(cr, 1, 0, 1); // Purple for text
cairo_text_extents_t extents;
cairo_text_extents(cr, cos_label, &extents);
double label_x = cos_x_position - extents.width / 2.0; // Center horizontally
double label_y = 0.05; // Adjust Y-coordinate to position below the line
cairo_move_to(cr, label_x, label_y);
cairo_show_text(cr, cos_label);
g_free(cos_label);
        
       
    }
}


if (draw_cos_wave2) {
    if (cos_x_position >= -1 && cos_x_position <= 1) {
        // Draw the cosine result dot
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the position
        cairo_arc(cr, cos_x_position, 0, 0.02, 0, 2 * M_PI); // Draw a small circle at the position
        cairo_fill(cr);

        // Draw a vertical line from the cosine result dot to the origin (0, 0)
        cairo_set_source_rgb(cr, 1, 0, 1); // Purple for the vertical line
        cairo_set_line_width(cr, 0.01); // Set line width to a thinner value
        cairo_move_to(cr, cos_x_position, 0); // Start from the dot position
        cairo_line_to(cr, 0, 0); // Draw the line to the origin (0, 0)
        cairo_stroke(cr);

        // Calculate the angle corresponding to the cosine value
        double angle = acos(cos_x_position); // Angle in radians (0 to π)
        double circle_x = cos(angle); // X-coordinate on the unit circle
        double circle_y = sin(angle); // Y-coordinate on the unit circle

        // Ensure the line is drawn in the upper part of the Cartesian plane
        if (circle_y < 0) {
            circle_y = -circle_y; // Reflect Y to make sure it is above the x-axis
        }

        // Draw a line from the point on the circle to the center
        cairo_set_source_rgb(cr, 0, 0, 1); // Blue for the circle-to-center line
        cairo_set_line_width(cr, 0.01); // Maintain the same line width
        cairo_move_to(cr, circle_x, circle_y); // Start from the calculated point on the circle
        cairo_line_to(cr, 0, 0); // Draw to the center
        cairo_stroke(cr);

        // Draw a vertical line upwards
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the vertical line
        cairo_set_line_width(cr, 0.01); // Maintain thinner line width
        cairo_move_to(cr, cos_x_position, 0); // Start from the dot position
        cairo_line_to(cr, circle_x, circle_y); // Stop at the circle point
        cairo_stroke(cr);
        
gchar *cos_label = g_strdup_printf("Cosine");
cairo_set_source_rgb(cr, 1, 0, 1); // Purple for text
cairo_text_extents_t extents;
cairo_text_extents(cr, cos_label, &extents);
double label_x = cos_x_position - extents.width / 2.0; // Center horizontally
double label_y = -0.05; // Adjust Y-coordinate to position below the line
cairo_move_to(cr, label_x, label_y);
cairo_show_text(cr, cos_label);
g_free(cos_label);        
    }
}

    // Draw the sine wave position
if (draw_sin_wave) {
    if (sin_y_position >= -1 && sin_y_position <= 1) {
        // Draw the sine result dot
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the position
        cairo_arc(cr, 0, sin_y_position, 0.02, 0, 2 * M_PI); // Draw a small circle at the position
        cairo_fill(cr);

        // Draw a vertical line from the sine result dot to the y-axis (origin)
        cairo_set_source_rgb(cr, 1, 0, 1); // Purple for the vertical line
        cairo_set_line_width(cr, 0.01); // Maintain thinner line width
        cairo_move_to(cr, 0, sin_y_position); // Start from the sine result on the y-axis
        cairo_line_to(cr, 0, 0); // Draw the line to the origin
        cairo_stroke(cr);

        // Calculate the angle corresponding to the sine value
        double angle = asin(sin_y_position); // Angle in radians (−π/2 to π/2)
        double circle_x = cos(angle); // X-coordinate on the unit circle
        double circle_y = sin(angle); // Y-coordinate on the unit circle

        // Draw the green line from the circle point to the sine dot on the y-axis
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the line
        cairo_set_line_width(cr, 0.01); // Maintain thinner line width
        cairo_move_to(cr, circle_x, circle_y); // Start at the circle point
        cairo_line_to(cr, 0, sin_y_position); // Stop at the sine result dot
        cairo_stroke(cr);

        // Draw a line from the point on the circle to the center
        cairo_set_source_rgb(cr, 0, 0, 1); // Blue for the circle-to-center line
        cairo_set_line_width(cr, 0.01); // Maintain the same line width
        cairo_move_to(cr, circle_x, circle_y); // Start from the calculated point on the circle
        cairo_line_to(cr, 0, 0); // Draw to the center
        cairo_stroke(cr);
        
gchar *sin_label = g_strdup_printf("Sine");
cairo_set_source_rgb(cr, 1, 0, 1); // Purple for text
cairo_move_to(cr, - 0.15, sin_y_position + 0.40);
cairo_show_text(cr, sin_label);
g_free(sin_label);        
    }
}



    // Draw the sine wave position
if (draw_sin_wave2) {
    if (sin_y_position >= -1 && sin_y_position <= 1) {
        // Draw the sine result dot
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the position
        cairo_arc(cr, 0, sin_y_position, 0.02, 0, 2 * M_PI); // Draw a small circle at the position
        cairo_fill(cr);

        // Draw a vertical line from the sine result dot to the y-axis (origin)
        cairo_set_source_rgb(cr, 1, 0, 1); // Purple for the vertical line
        cairo_set_line_width(cr, 0.01); // Maintain thinner line width
        cairo_move_to(cr, 0, sin_y_position); // Start from the sine result on the y-axis
        cairo_line_to(cr, 0, 0); // Draw the line to the origin
        cairo_stroke(cr);

        // Calculate the angle corresponding to the sine value
        double angle = asin(sin_y_position); // Angle in radians (−π/2 to π/2)
        angle = M_PI - angle; // Adjust angle for the second quadrant (90°–180°)

        double circle_x = cos(angle); // X-coordinate on the unit circle
        double circle_y = sin(angle); // Y-coordinate on the unit circle

        // Draw the green line from the circle point to the sine dot on the y-axis
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the line
        cairo_set_line_width(cr, 0.01); // Maintain thinner line width
        cairo_move_to(cr, circle_x, circle_y); // Start at the circle point
        cairo_line_to(cr, 0, sin_y_position); // Stop at the sine result dot
        cairo_stroke(cr);

        // Draw a line from the point on the circle to the center
        cairo_set_source_rgb(cr, 0, 0, 1); // Blue for the circle-to-center line
        cairo_set_line_width(cr, 0.01); // Maintain the same line width
        cairo_move_to(cr, circle_x, circle_y); // Start from the calculated point on the circle
        cairo_line_to(cr, 0, 0); // Draw to the center
        cairo_stroke(cr);

gchar *sin_label = g_strdup_printf("Sine");
cairo_set_source_rgb(cr, 1, 0, 1); // Purple for text
cairo_move_to(cr, + 0.15, sin_y_position + 0.40);
cairo_show_text(cr, sin_label);
g_free(sin_label); 

    }
}


   // Draw the sine wave position
if (draw_sin_wave3) {
    if (sin_y_position >= -1 && sin_y_position <= 1) {
        // Draw the sine result dot
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the position
        cairo_arc(cr, 0, sin_y_position, 0.02, 0, 2 * M_PI); // Draw a small circle at the position
        cairo_fill(cr);

        // Draw a vertical line from the sine result dot to the y-axis (origin)
        cairo_set_source_rgb(cr, 1, 0, 1); // Purple for the vertical line
        cairo_set_line_width(cr, 0.01); // Maintain thinner line width
        cairo_move_to(cr, 0, sin_y_position); // Start from the sine result on the y-axis
        cairo_line_to(cr, 0, 0); // Draw the line to the origin
        cairo_stroke(cr);

        // Calculate the angle corresponding to the sine value
        double angle = asin(sin_y_position); // Angle in radians (−π/2 to π/2)
        angle = M_PI - angle; // Shift the angle to the third quadrant (180°–270°)

        // Compute the circle coordinates
        double circle_x = cos(angle); // Negative x-coordinate for the third quadrant
        double circle_y = sin(angle); // Negative y-coordinate for the third quadrant

        // Draw the green line from the circle point to the sine dot on the y-axis
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the line
        cairo_set_line_width(cr, 0.01); // Maintain thinner line width
        cairo_move_to(cr, circle_x, circle_y); // Start at the circle point
        cairo_line_to(cr, 0, sin_y_position); // Stop at the sine result dot
        cairo_stroke(cr);

        // Draw a line from the point on the circle to the center
        cairo_set_source_rgb(cr, 0, 0, 1); // Blue for the circle-to-center line
        cairo_set_line_width(cr, 0.01); // Maintain the same line width
        cairo_move_to(cr, circle_x, circle_y); // Start from the calculated point on the circle
        cairo_line_to(cr, 0, 0); // Draw to the center
        cairo_stroke(cr);
        
gchar *sin_label = g_strdup_printf("Sine");
cairo_set_source_rgb(cr, 1, 0, 1); // Purple for text
cairo_move_to(cr, + 0.15, sin_y_position - 0.40);
cairo_show_text(cr, sin_label);
g_free(sin_label);         
    }
}

   // Draw the sine wave position
if (draw_sin_wave4) {
    if (sin_y_position >= -1 && sin_y_position <= 1) {
        // Draw the sine result dot
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the position
        cairo_arc(cr, 0, sin_y_position, 0.02, 0, 2 * M_PI); // Draw a small circle at the position
        cairo_fill(cr);

        // Draw a vertical line from the sine result dot to the y-axis (origin)
        cairo_set_source_rgb(cr, 1, 0, 1); // Purple for the vertical line
        cairo_set_line_width(cr, 0.01); // Maintain thinner line width
        cairo_move_to(cr, 0, sin_y_position); // Start from the sine result on the y-axis
        cairo_line_to(cr, 0, 0); // Draw the line to the origin
        cairo_stroke(cr);

        // Calculate the angle corresponding to the sine value
        double angle = asin(sin_y_position); // Angle in radians (−π/2 to π/2)
        angle = 2 * M_PI + angle; // Shift the angle to the fourth quadrant (270°–360°)

        // Compute the circle coordinates
        double circle_x = cos(angle); // Positive x-coordinate for the fourth quadrant
        double circle_y = sin(angle); // Negative y-coordinate for the fourth quadrant

        // Draw the green line from the circle point to the sine dot on the y-axis
        cairo_set_source_rgb(cr, 0, 1, 0); // Green for the line
        cairo_set_line_width(cr, 0.01); // Maintain thinner line width
        cairo_move_to(cr, circle_x, circle_y); // Start at the circle point
        cairo_line_to(cr, 0, sin_y_position); // Stop at the sine result dot
        cairo_stroke(cr);

        // Draw a line from the point on the circle to the center
        cairo_set_source_rgb(cr, 0, 0, 1); // Blue for the circle-to-center line
        cairo_set_line_width(cr, 0.01); // Maintain the same line width
        cairo_move_to(cr, circle_x, circle_y); // Start from the calculated point on the circle
        cairo_line_to(cr, 0, 0); // Draw to the center
        cairo_stroke(cr);
 
gchar *sin_label = g_strdup_printf("Sine");
cairo_set_source_rgb(cr, 1, 0, 1); // Purple for text
cairo_move_to(cr, - 0.15, sin_y_position - 0.40);
cairo_show_text(cr, sin_label);
g_free(sin_label);  
        
    }
}

    cairo_restore(cr); // Restore the saved state
}
    
if ((draw_tan || draw_atan)) {

 cairo_save(cr); 
 
        for (double i = -1.0; i <= 1.0; i += 0.5) {
            char label[10];
            snprintf(label, sizeof(label), "%.1f", i);
            cairo_move_to(cr, i - 0.05, 0.1);
            cairo_show_text(cr, label);
            cairo_move_to(cr, 0.05, -i + 0.05);
            cairo_show_text(cr, label);

  double calculated_tan = expression_result;
  
//DRAW TAN  
if (draw_tan) {
cairo_save(cr);
cairo_set_operator(cr, CAIRO_OPERATOR_CLEAR);
cairo_paint(cr);
cairo_restore(cr);


    cairo_set_source_rgba(cr, 0, 1, 0, 1.0); // Green for the position with full opacity

    // Calculate the coordinates for the circle
    double x_coord = width / 2 + origin * 20.0; // Scaling origin by 20 (grid step size)
    double y_coord = height / 2 - calculated_tan * 20.0; // Scaling y by 20
    double radius = 5.0; // Radius of the circle
   
    // Draw a small circle at the position
    cairo_arc(cr, x_coord, y_coord, radius, 0, 2 * M_PI);
    cairo_fill(cr);
    
    cairo_new_path(cr);   
    cairo_set_source_rgba(cr, 1, 0, 1, 1.0); // Blue for the tan curve with full opacity
    cairo_set_line_width(cr, 1);

    bool is_first_point = true;

    // Loop to draw the tangent curve only in the range -12 to 12
   
    for (double x = -1.6; x <= 1.6; x += 0.01) {
        double y = tan(x);
        if (fabs(y) < 10) { // This avoids extremely large values that would distort the curve
            if (is_first_point) {
                cairo_move_to(cr, width / 2 + x * 20, height / 2 - y * 20); // Adjust x-coordinate scaling
                is_first_point = false;
            } else {
                cairo_line_to(cr, width / 2 + x * 20, height / 2 - y * 20); // Adjust x-coordinate scaling
            }
        } else {
            is_first_point = true; // Reset for the next valid segment
        }
    }
    cairo_stroke(cr);
    cairo_restore(cr); // Restore to the previous state
}
//END TAN
//START ATAN
    if (draw_atan) {

    cairo_save(cr);
    cairo_set_operator(cr, CAIRO_OPERATOR_CLEAR);
    cairo_paint(cr);
    cairo_restore(cr);

    cairo_set_source_rgba(cr, 0, 1, 0, 1.0); // Green for the position with full opacity

    // Calculate the coordinates for the circle
    double x_coord = width / 2 + origin * 20.0; // Scaling origin by 20 (grid step size)
    double y_coord = height / 2 - calculated_atan * 20.0; // Scaling y by 20
    double radius = 5.0; // Radius of the circle
   
    // Draw a small circle at the position
    cairo_arc(cr, x_coord, y_coord, radius, 0, 2 * M_PI);
    cairo_fill(cr);
    
    cairo_new_path(cr);   
    cairo_set_source_rgba(cr, 1, 0, 1, 1.0); // Blue for the tan curve with full opacity
    cairo_set_line_width(cr, 1);

    bool is_first_point = true;
    cairo_move_to(cr, width / 2 + (-12) * 20.0, height / 2 - atan(-12) * 20.0); // Adjust initial move_to point

    for (double x = -12; x <= 12; x += 0.01) { // Adjust loop range and increment
        double y = atan(x);
        cairo_line_to(cr, width / 2 + x * 20.0, height / 2 - y * 20.0); // Scaling x and y by 20
            }
            cairo_stroke(cr);
            cairo_restore(cr); // Restore to the previous state
        }
}
//END ATAN

}}

// Declare the drawing area globally
GtkWidget *calcdrawing_area;

// Function to set up the drawing area
void setup_drawing_area(GtkWidget *calculatorbox) {
    calcdrawing_area = gtk_drawing_area_new();
    gtk_widget_set_size_request(calcdrawing_area, 400, 400);
    gtk_widget_set_name(calcdrawing_area, "calcdrawing_area");
    g_signal_connect(G_OBJECT(calcdrawing_area), "draw", G_CALLBACK(draw_callback), NULL);
    gtk_box_pack_start(GTK_BOX(calculatorbox), calcdrawing_area, TRUE, TRUE, 0);
}

// Function to update the entry expression
static void on_number_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *button_label = gtk_button_get_label(GTK_BUTTON(button));
    
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    gchar *new_text = g_strdup_printf("%s%s", current_text, button_label);
    gtk_entry_set_text(entry_expression, new_text);
    g_free(new_text);
}

// operation +-*/
static void on_plus_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);

    // Use previous result if available, otherwise take input number
    if (expression_result != 0) {
        first_number = expression_result;
    } else {
        first_number = atof(current_text);
    }

    // Clear input for new number entry
    gtk_entry_set_text(entry_expression, "");
    operation = 1; // Set operation to addition
}

static void on_subtraction_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);

    // Use previous result if available, otherwise take input number
    if (expression_result != 0) {
        first_number = expression_result;
    } else {
        first_number = atof(current_text);
    }

    // Clear input for new number entry
    gtk_entry_set_text(entry_expression, "");
    operation = 2; // Set operation to subtraction
}

static void on_multiplication_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);

    // If an operation was previously performed, use the last result
    if (expression_result != 0) {
        first_number = expression_result;
    } else {
        first_number = atof(current_text);
    }

    gtk_entry_set_text(entry_expression, ""); // Clear input for new number
    operation = 3; // Set operation to multiplication
}


static void on_division_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    first_number = atof(current_text);
    gtk_entry_set_text(entry_expression, "");
    operation = 4; // Set operation to division
}

//end operation +-*/

//XROOT
static void on_xroot_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    first_number = atof(current_text);
    gtk_entry_set_text(entry_expression, "");
    operation = 5; 
}

//XPOWER
static void on_xpower_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    first_number = atof(current_text);
    gtk_entry_set_text(entry_expression, "");
    operation = 6; 
}

// EQUAL BUTTON
static void on_equal_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);

    double second_number = atof(current_text);

    // Use previous result if no new input is entered
    if (strlen(current_text) == 0) {
        second_number = expression_result;
    }

    if (operation == 1) { // Addition
        expression_result = first_number + second_number;
    } else if (operation == 2) { // Subtraction
        expression_result = first_number - second_number;
    } else if (operation == 3) { // Multiplication
        expression_result = first_number * second_number;
    } else if (operation == 4) { // Division
        if (second_number != 0) {
            expression_result = first_number / second_number;
        } else {
            // Handle division by zero
            GtkEntry *entry_result = GTK_ENTRY(g_object_get_data(G_OBJECT(button), "entry_result"));
            gtk_entry_set_text(entry_result, "Error");
            return;
        }
    } else if (operation == 5) { // x-th root
        if (second_number == 0) {
            // Handle root of zero
            GtkEntry *entry_result = GTK_ENTRY(g_object_get_data(G_OBJECT(button), "entry_result"));
            gtk_entry_set_text(entry_result, "Error: Root cannot be zero");
            return;
        } else {
            expression_result = pow(first_number, 1.0 / second_number);
        }
    } else if (operation == 6) { // Exponentiation
        if (second_number == 0) {
            expression_result = 1;
        } else {
            expression_result = pow(first_number, second_number);
        }
    }

    // Display the result
    GtkEntry *entry_result = GTK_ENTRY(g_object_get_data(G_OBJECT(button), "entry_result"));
    gchar result_text[32];
    snprintf(result_text, sizeof(result_text), "%.2f", expression_result);
    gtk_entry_set_text(entry_result, result_text);

    // Update first_number to store the result for continued calculations
    first_number = expression_result;

    // Don't reset operation immediately; keep it active for further calculations
}


// EULERS
static void on_eulero_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    gdouble number = atof(current_text);
    expression_result = exp(number);  // Update expression_result
    gchar result_text[64];
    g_snprintf(result_text, sizeof(result_text), "%.6f", expression_result);
    gtk_entry_set_text(entry_expression, result_text);
    
    // Queue the draw after updating expression_result
    gtk_widget_queue_draw(calcdrawing_area);

}


// NEGATIVE POWER
static void on_negativepower_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);

    gdouble number = atof(current_text);
    gdouble result = (number != 0) ? 1.0 / number : 0.0; // Avoid division by zero

    gchar result_text[64]; 
    g_snprintf(result_text, sizeof(result_text), "%.6f", result);
    gtk_entry_set_text(entry_expression, result_text);
    
    // Queue the draw after updating expression_result
    expression_result = result;  // Update expression_result with the result

    gtk_widget_queue_draw(calcdrawing_area);

}


//SQUAREROOT
static void on_squareroot_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    double number = atof(current_text);
    
    if (number < 0) {
        GtkEntry *entry_result = GTK_ENTRY(g_object_get_data(G_OBJECT(button), "entry_result"));
        gtk_entry_set_text(entry_result, "Error");
    } else {
        double result = sqrt(number);

        GtkEntry *entry_result = GTK_ENTRY(g_object_get_data(G_OBJECT(button), "entry_result"));
        gchar result_text[32];
        snprintf(result_text, sizeof(result_text), "%.2f", result);
        gtk_entry_set_text(entry_result, result_text);
        expression_result = result;  // Update expression_result with the result
    }  
    // Queue the draw after updating expression_result
    gtk_widget_queue_draw(calcdrawing_area);
}

//CUBEROOT
static void on_cuberoot_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);

    double number = atof(current_text);
    
    if (number < 0) {
        GtkEntry *entry_result = GTK_ENTRY(g_object_get_data(G_OBJECT(button), "entry_result"));
        gtk_entry_set_text(entry_result, "Error");
    } else {
        double result = cbrt(number);
        
        GtkEntry *entry_result = GTK_ENTRY(g_object_get_data(G_OBJECT(button), "entry_result"));
        gchar result_text[32];
        snprintf(result_text, sizeof(result_text), "%.2f", result);
        gtk_entry_set_text(entry_result, result_text);
        expression_result = result;  // Update expression_result with the result
        // Queue the draw after updating expression_result
        gtk_widget_queue_draw(calcdrawing_area);
    }

}

//SQUARE POWER
static void on_squarepower_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    gdouble number = atof(current_text);
    gdouble result = number * number;
    gchar result_text[64];
    g_snprintf(result_text, sizeof(result_text), "%.6f", result);
    gtk_entry_set_text(entry_expression, result_text);
    // Queue the draw after updating expression_result
    expression_result = result;  // Update expression_result with the result
    gtk_widget_queue_draw(calcdrawing_area);
}

//CUBE POWER
static void on_cubepower_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    gdouble number = atof(current_text);
    gdouble result = number * number * number;
    gchar result_text[64]; 
    g_snprintf(result_text, sizeof(result_text), "%.6f", result); 
    gtk_entry_set_text(entry_expression, result_text);
    // Queue the draw after updating expression_result
    expression_result = result;  // Update expression_result with the result
    gtk_widget_queue_draw(calcdrawing_area);
}

//PGRECO
static void on_pgreco_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    gdouble number = atof(current_text);
    expression_result = number * 3.141592;
    gchar result_text[64]; 
    g_snprintf(result_text, sizeof(result_text), "%.6f", expression_result);
    gtk_entry_set_text(entry_expression, result_text);
    
  //  g_print("Result calculated: %.6f\n", expression_result);
    
    // Queue the draw after updating expression_result
    gtk_widget_queue_draw(calcdrawing_area);

}


//LOG
static void on_log_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    gdouble number = atof(current_text);
    gdouble result = log(number);
    gchar result_text[64]; 
    g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

    gtk_entry_set_text(entry_expression, result_text);
    // Queue the draw after updating expression_result
    expression_result = result;  // Update expression_result with the result
    gtk_widget_queue_draw(calcdrawing_area);
}


// Converts degrees to radians
static double degrees_to_radians(double degrees) {
    return degrees * (M_PI / 180.0);
}


// Declare the variables at a global scope if they are used across different functions
//COS 
static void on_cos_button_clicked(GtkWidget *button, gpointer user_data) { 
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);

    gdouble number = atof(current_text);
    if (number >= 0 && number <= 180) { // Check if the number is between 0 and 180
        gdouble radians = degrees_to_radians(number); // Convert degrees to radians
        gdouble result = cos(radians);
        gchar result_text[64]; 
        g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

        gtk_entry_set_text(entry_expression, result_text); 
        // Queue the draw after updating expression_result 
        expression_result = result; // Update expression_result with the result 
        calculated_cosine = result; // Update the calculated_cosine with the result
        draw_cos_wave = TRUE; // Set the flag to draw the cosine wave 
        gtk_widget_queue_draw(calcdrawing_area); 
    } else if (number > 180 && number <= 360) { // Check if the number is between 180 and 360
        gdouble radians = degrees_to_radians(number); // Convert degrees to radians
        gdouble result = cos(radians);
        gchar result_text[64]; 
        g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

        gtk_entry_set_text(entry_expression, result_text); 
        // Queue the draw after updating expression_result 
        expression_result = result; // Update expression_result with the result 
        calculated_cosine = result; // Update the calculated_cosine with the result
        draw_cos_wave2 = TRUE; // Set the flag to draw the cosine wave 
        gtk_widget_queue_draw(calcdrawing_area); 
    } else {
        // Optionally, handle cases where the number is out of the desired range
        gtk_entry_set_text(entry_expression, "Number out of range (0-360)");
    }
}

//ACOS
static void on_acos_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    gdouble number = atof(current_text);
    gdouble result = acos(number);
    gchar result_text[64]; 
    g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

    gtk_entry_set_text(entry_expression, result_text);
    // Queue the draw after updating expression_result
    expression_result = result;  // Update expression_result with the result
    gtk_widget_queue_draw(calcdrawing_area);
}

//COSH
static void on_cosh_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    gdouble number = atof(current_text);
    gdouble result = cosh(number);
    gchar result_text[64]; 
    g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

    gtk_entry_set_text(entry_expression, result_text);
    // Queue the draw after updating expression_result
    expression_result = result;  // Update expression_result with the result
    gtk_widget_queue_draw(calcdrawing_area);
}

//SIN
static void on_sin_button_clicked(GtkWidget *button, gpointer user_data) { 
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);

    gdouble number = atof(current_text);
    if (number >= 0 && number <= 90) { // Check if the number is between 0 and 180
        gdouble radians = degrees_to_radians(number); // Convert degrees to radians
        gdouble result = sin(radians);
        gchar result_text[64]; 
        g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

        gtk_entry_set_text(entry_expression, result_text); 
        // Queue the draw after updating expression_result 
        expression_result = result; // Update expression_result with the result 
        calculated_sin = result; // Update the calculated_cosine with the result
        draw_sin_wave = TRUE; // Set the flag to draw the cosine wave 
        gtk_widget_queue_draw(calcdrawing_area); 
    }  else if (number > 90 && number <= 180) { // Check if the number is between 180 and 360
        gdouble radians = degrees_to_radians(number); // Convert degrees to radians
        gdouble result = sin(radians);
        gchar result_text[64]; 
        g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

        gtk_entry_set_text(entry_expression, result_text); 
        // Queue the draw after updating expression_result 
        expression_result = result; // Update expression_result with the result 
        calculated_sin = result; // Update the calculated_cosine with the result
        draw_sin_wave2 = TRUE; // Set the flag to draw the cosine wave 
        gtk_widget_queue_draw(calcdrawing_area); 
    } 
    else if (number > 180 && number <= 270) { // Check if the number is between 180 and 360
        gdouble radians = degrees_to_radians(number); // Convert degrees to radians
        gdouble result = sin(radians);
        gchar result_text[64]; 
        g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

        gtk_entry_set_text(entry_expression, result_text); 
        // Queue the draw after updating expression_result 
        expression_result = result; // Update expression_result with the result 
        calculated_sin = result; // Update the calculated_cosine with the result
        draw_sin_wave3 = TRUE; // Set the flag to draw the cosine wave 
        gtk_widget_queue_draw(calcdrawing_area); 
    }  
    else if (number > 270 && number <= 360) { // Check if the number is between 180 and 360
        gdouble radians = degrees_to_radians(number); // Convert degrees to radians
        gdouble result = sin(radians);
        gchar result_text[64]; 
        g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

        gtk_entry_set_text(entry_expression, result_text); 
        // Queue the draw after updating expression_result 
        expression_result = result; // Update expression_result with the result 
        calculated_sin = result; // Update the calculated_cosine with the result
        draw_sin_wave4 = TRUE; // Set the flag to draw the cosine wave 
        gtk_widget_queue_draw(calcdrawing_area); 
    } else {
        // Optionally, handle cases where the number is out of the desired range
        gtk_entry_set_text(entry_expression, "Number out of range (0-360)");
    }
}

//ASIN
static void on_asin_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    gdouble number = atof(current_text);
    gdouble result = asin(number);
    gchar result_text[64]; 
    g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

    gtk_entry_set_text(entry_expression, result_text);
    // Queue the draw after updating expression_result
    expression_result = result;  // Update expression_result with the result

    gtk_widget_queue_draw(calcdrawing_area);
}
//SINH
static void on_sinh_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    gdouble number = atof(current_text);
    gdouble result = sinh(number);
    gchar result_text[64]; 
    g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

    gtk_entry_set_text(entry_expression, result_text);
    // Queue the draw after updating expression_result
    expression_result = result;  // Update expression_result with the result
    gtk_widget_queue_draw(calcdrawing_area);
}
//TAN
static void on_tan_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);

    gdouble number = atof(current_text); // Get the input number
    gdouble result = tan(number); // Calculate the tangent
    gchar result_text[64];
    g_snprintf(result_text, sizeof(result_text), "%.6f", result); // Format the result

    gtk_entry_set_text(entry_expression, result_text); // Update the entry with the result
    expression_result = result;  // Store the result
    calculated_tan = result; // Set the calculated_tan variable
    origin = number;
    draw_tan = TRUE; // Enable drawing the tangent curve
    gtk_widget_queue_draw(calcdrawing_area); // Queue the draw
}

//ATAN
static void on_atan_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    const gchar *current_text = gtk_entry_get_text(entry_expression);
    
    gdouble number = atof(current_text);
    gdouble result = atan(number);
    gchar result_text[64]; 
    g_snprintf(result_text, sizeof(result_text), "%.6f", result); 

    gtk_entry_set_text(entry_expression, result_text);
    // Queue the draw after updating expression_result
    expression_result = result;  // Update expression_result with the result
    calculated_atan = result;
    origin = number;
    draw_atan = TRUE; // Set the flag to draw the cosine wave 
    gtk_widget_queue_draw(calcdrawing_area); 
}

// AC
static void on_ac_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkEntry *entry_expression = GTK_ENTRY(user_data);
    GtkEntry *entry_result = GTK_ENTRY(g_object_get_data(G_OBJECT(button), "entry_result"));
    
    // Clear the entries
    gtk_entry_set_text(entry_expression, "");
    gtk_entry_set_text(entry_result, "");
    
    // Reset the expression result
    expression_result = 0.0;
  
    // Reset the cosine wave flag
    draw_cos_wave = FALSE;
    draw_cos_wave2 = FALSE;   
    draw_sin_wave = FALSE;
    draw_sin_wave2 = FALSE;  
    draw_sin_wave3 = FALSE; 
    draw_sin_wave4 = FALSE; 
    draw_atan = FALSE; 
    draw_tan = FALSE;
    draw_parabolay = FALSE;
    draw_parabolax = FALSE;
    draw_line = FALSE;
    draw_eq2 = FALSE;
    draw_eq3 = FALSE;
    draw_archim = FALSE;
    draw_expo = FALSE;
    draw_proportions = FALSE;
    draw_ineq2 = FALSE;
    // Optionally, redraw the drawing area if needed
    GtkWidget *calcdrawing_area = GTK_WIDGET(g_object_get_data(G_OBJECT(button), "calcdrawing_area"));
    gtk_widget_queue_draw(calcdrawing_area);
}


void on_parabola_button_clicked(GtkWidget *button, GtkWidget *parabolabox) {
    if (gtk_widget_get_visible(parabolabox)) {
        gtk_widget_hide(parabolabox);
    } else {
        gtk_widget_show(parabolabox);
    }
}

void on_show_parabolay_button_clicked(GtkButton *button, gpointer user_data) {
    static gboolean is_empty = FALSE;
    gpointer *data = (gpointer *)user_data;
    GtkWidget **entries = (GtkWidget **)data[0];
    GtkWidget *label_result_parabolay = GTK_WIDGET(data[1]);
    GtkWidget *label_vertex_x_parabolay = GTK_WIDGET(data[2]);
    GtkWidget *label_delta_parabolay = GTK_WIDGET(data[3]);
    GtkWidget *label_vertex_y_parabolay = GTK_WIDGET(data[4]);
    GtkWidget *label_focus_x_parabolay = GTK_WIDGET(data[5]);
    GtkWidget *label_focus_y_parabolay = GTK_WIDGET(data[6]);
    GtkWidget *sim_axys_parabolay = GTK_WIDGET(data[7]);
    GtkWidget *directrix = GTK_WIDGET(data[8]);

    if (is_empty) {
        // Reset all values
        for (int i = 0; i < 3; i++) {
            gtk_entry_set_text(GTK_ENTRY(entries[i]), "");
        }
        gtk_label_set_text(GTK_LABEL(label_result_parabolay), "");
        gtk_label_set_text(GTK_LABEL(label_vertex_x1_parabolay), "");
        gtk_label_set_text(GTK_LABEL(label_vertex_x2_parabolay), "");
        gtk_label_set_text(GTK_LABEL(label_delta_parabolay), "");
        gtk_label_set_text(GTK_LABEL(label_vertex_y_parabolay), "");
        gtk_label_set_text(GTK_LABEL(label_focus_x_parabolay), "");
        gtk_label_set_text(GTK_LABEL(label_focus_y_parabolay), "");
        gtk_label_set_text(GTK_LABEL(sim_axys_parabolay), "");
        gtk_label_set_text(GTK_LABEL(directrix), "");
        gtk_button_set_label(button, "Show Parabola");
    } else {
        // Check if entries are valid
        for (int i = 0; i < 3; i++) {
            if (!GTK_IS_ENTRY(entries[i])) {
                //g_print("Invalid GtkEntry pointer at index %d\n", i);
                return;
            }
        }

        const gchar *a_value = gtk_entry_get_text(GTK_ENTRY(entries[0]));
        const gchar *b_value = gtk_entry_get_text(GTK_ENTRY(entries[1]));
        const gchar *c_value = gtk_entry_get_text(GTK_ENTRY(entries[2]));

        gchar *result_text = g_strdup_printf("y = %sx² + %sx + %s", a_value, b_value, c_value);
        gtk_label_set_text(GTK_LABEL(label_result_parabolay), result_text);
        g_free(result_text);

        double a = atof(a_value);
        double b = atof(b_value);
        double c = atof(c_value);

        double delta = (b * b) - (4 * a * c);
        gchar *delta_text = g_strdup_printf("Δ = %.2f", delta);
        gtk_label_set_text(GTK_LABEL(label_delta_parabolay), delta_text);
        g_free(delta_text);
        py_delta = delta; 
        double x1_vertex = (-b + sqrt(delta)) / (2 * a);
        gchar *vertex_text = g_strdup_printf("x1 vertex = %.2f", x1_vertex);
        gtk_label_set_text(GTK_LABEL(label_vertex_x1_parabolay), vertex_text);
        g_free(vertex_text);
        py_xvertex = x1_vertex; 
        
        double x2_vertex = (-b - sqrt(delta)) / (2 * a);
        gchar *vertex2_text = g_strdup_printf("x2 vertex = %.2f", x2_vertex);
        gtk_label_set_text(GTK_LABEL(label_vertex_x2_parabolay), vertex2_text);
        g_free(vertex2_text);
        py_x2vertex = x2_vertex; 
        
        double y_vertex = -delta / (4 * a);
        gchar *yvertex_text = g_strdup_printf("y vertex = %.2f", y_vertex);
        gtk_label_set_text(GTK_LABEL(label_vertex_y_parabolay), yvertex_text);
        g_free(yvertex_text);
        py_yvertex = y_vertex;
        double x_focus = -b / (2 * a);
        gchar *focus_text = g_strdup_printf("x focus = %.2f", x_focus);
        gtk_label_set_text(GTK_LABEL(label_focus_x_parabolay), focus_text);
        g_free(focus_text);
        py_xfocus = x_focus;
        double y_focus = 1 - delta / (4 * a);
        gchar *focusy_text = g_strdup_printf("y focus = %.2f", y_focus);
        gtk_label_set_text(GTK_LABEL(label_focus_y_parabolay), focusy_text);
        g_free(focusy_text);
        py_yfocus = y_focus;
        double sim_axys = -b / (2 * a);
        gchar *sim_axys_text = g_strdup_printf("Axis of Symmetry x= %.2f", sim_axys);
        gtk_label_set_text(GTK_LABEL(sim_axys_parabolay), sim_axys_text);
        g_free(sim_axys_text);
        py_simaxys = sim_axys;
        double dir = 1 + delta / (4 * a);
        gchar *dir_text = g_strdup_printf("Directrix y= %.2f", dir);
        gtk_label_set_text(GTK_LABEL(directrix), dir_text);
        g_free(dir_text);
        py_dir = dir;
        gtk_button_set_label(button, "Clear");
        
        draw_parabolay = TRUE; // Set flag to true
       // g_print("Debug: draw_parabolay set to TRUE\n"); // Debug statement

        if (GTK_IS_WIDGET(calcdrawing_area)) {
            gtk_widget_queue_draw(calcdrawing_area); // Force redraw of the drawing area
        } else {
            //g_print("Error: drawing_area is not a valid GtkWidget\n");
        }
    }


    is_empty = !is_empty;
}


///PARABOLA X
void on_parabolax_button_clicked(GtkWidget *button, GtkWidget *parabola_x_box) {
    if (gtk_widget_get_visible(parabola_x_box)) {
        gtk_widget_hide(parabola_x_box);
    } else {
        gtk_widget_show(parabola_x_box);
    }
}

void on_show_parabolax_button_clicked(GtkButton *button, gpointer user_data) {
    static gboolean is_empty = FALSE;
    gpointer *datax = (gpointer *)user_data;
    GtkWidget **entriesx = (GtkWidget **)datax[0];
    GtkWidget *label_result_parabolax = GTK_WIDGET(datax[1]);
    GtkWidget *label_vertex_x_parabolax = GTK_WIDGET(datax[2]);
    GtkWidget *label_delta_parabolax = GTK_WIDGET(datax[3]);
    GtkWidget *label_vertex_y1_parabolax = GTK_WIDGET(datax[4]);
    GtkWidget *label_vertex_y2_parabolax = GTK_WIDGET(datax[5]);
    GtkWidget *label_focus_x_parabolax = GTK_WIDGET(datax[6]);
    GtkWidget *label_focus_y_parabolax = GTK_WIDGET(datax[7]);
    GtkWidget *sim_axys_parabolax = GTK_WIDGET(datax[8]);
    GtkWidget *directrix_x = GTK_WIDGET(datax[9]);

    if (is_empty) {
        // Reset all values
        for (int i = 0; i < 3; i++) {
            gtk_entry_set_text(GTK_ENTRY(entriesx[i]), "");
        }
        gtk_label_set_text(GTK_LABEL(label_result_parabolax), "");
        gtk_label_set_text(GTK_LABEL(label_vertex_x_parabolax), "");
        gtk_label_set_text(GTK_LABEL(label_delta_parabolax), "");
        gtk_label_set_text(GTK_LABEL(label_vertex_y1_parabolax), "");
        gtk_label_set_text(GTK_LABEL(label_focus_x_parabolax), "");
        gtk_label_set_text(GTK_LABEL(label_focus_y_parabolax), "");
        gtk_label_set_text(GTK_LABEL(sim_axys_parabolax), "");
        gtk_label_set_text(GTK_LABEL(directrix_x), "");
        gtk_button_set_label(button, "Show Parabola");
    } else {
        // Check if entries are valid
        for (int i = 0; i < 3; i++) {
            if (!GTK_IS_ENTRY(entriesx[i])) {
                //g_print("Invalid GtkEntry pointer at index %d\n", i);
                return;
            }
        }

        const gchar *a_value = gtk_entry_get_text(GTK_ENTRY(entriesx[0]));
        const gchar *b_value = gtk_entry_get_text(GTK_ENTRY(entriesx[1]));
        const gchar *c_value = gtk_entry_get_text(GTK_ENTRY(entriesx[2]));

        gchar *result_text = g_strdup_printf("x = %sy² + %sy + %s", a_value, b_value, c_value);
        gtk_label_set_text(GTK_LABEL(label_result_parabolax), result_text);
        g_free(result_text);

        double a = atof(a_value);
        double b = atof(b_value);
        double c = atof(c_value);

        double delta = (b * b) - (4 * a * c);
        gchar *delta_text = g_strdup_printf("Δ = %.2f", delta);
        gtk_label_set_text(GTK_LABEL(label_delta_parabolax), delta_text);
        g_free(delta_text);
        py_delta = delta; 
        double x_vertex = -delta / (4 * a);
        gchar *vertex_text = g_strdup_printf("x vertex = %.2f", x_vertex);
        gtk_label_set_text(GTK_LABEL(label_vertex_x_parabolax), vertex_text);
        g_free(vertex_text);
        py_xvertex = x_vertex; 
        double y1_vertex = (-b + sqrt(delta)) / (2 * a);
        gchar *yvertex1_text = g_strdup_printf("y1 vertex = %.2f", y1_vertex);
        gtk_label_set_text(GTK_LABEL(label_vertex_y1_parabolax), yvertex1_text);
        g_free(yvertex1_text);
        py_yvertex = y1_vertex;
        
        double y2_vertex = (-b - sqrt(delta)) / (2 * a);
        gchar *yvertex2_text = g_strdup_printf("y2 vertex = %.2f", y2_vertex);
        gtk_label_set_text(GTK_LABEL(label_vertex_y2_parabolax), yvertex2_text);
        g_free(yvertex2_text);
        py_y2vertex = y2_vertex;        
        
        double x_focus = -b / (2 * a);
        gchar *focus_text = g_strdup_printf("x focus = %.2f", x_focus);
        gtk_label_set_text(GTK_LABEL(label_focus_x_parabolax), focus_text);
        g_free(focus_text);
        py_xfocus = x_focus;
        double y_focus = 1 - delta / (4 * a);
        gchar *focusy_text = g_strdup_printf("y focus = %.2f", y_focus);
        gtk_label_set_text(GTK_LABEL(label_focus_y_parabolax), focusy_text);
        g_free(focusy_text);
        py_yfocus = y_focus;
        double sim_axys = -b / (2 * a);
        gchar *sim_axys_text = g_strdup_printf("Axis of Symmetry x= %.2f", sim_axys);
        gtk_label_set_text(GTK_LABEL(sim_axys_parabolax), sim_axys_text);
        g_free(sim_axys_text);
        py_simaxys = sim_axys;
        double dir = 1 + delta / (4 * a);
        gchar *dir_text = g_strdup_printf("Directrix y= %.2f", dir);
        gtk_label_set_text(GTK_LABEL(directrix_x), dir_text);
        g_free(dir_text);
        py_dir = dir;
        gtk_button_set_label(button, "Clear");
        
        draw_parabolax = TRUE; // Set flag to true
       // g_print("Debug: draw_parabolay set to TRUE\n"); // Debug statement

        if (GTK_IS_WIDGET(calcdrawing_area)) {
            gtk_widget_queue_draw(calcdrawing_area); // Force redraw of the drawing area
        } else {
            //g_print("Error: drawing_area is not a valid GtkWidget\n");
        }
    }


    is_empty = !is_empty;
}

//LINEAR EQUATION
void on_line_button_clicked(GtkWidget *button, GtkWidget *line_box) {
    if (gtk_widget_get_visible(line_box)) {
        gtk_widget_hide(line_box);
    } else {
        gtk_widget_show(line_box);
    }
}

void on_show_line_button_clicked(GtkButton *button, gpointer user_data) {
    static gboolean is_empty = FALSE;
    gpointer *data_line = (gpointer *)user_data;
    GtkWidget **entries_line = (GtkWidget **)data_line[0];
    GtkWidget *label_m_line = GTK_WIDGET(data_line[1]);
    GtkWidget *label_b_line = GTK_WIDGET(data_line[2]);
    GtkWidget *label_result_line = GTK_WIDGET(data_line[3]);

    if (is_empty) {
        // Reset all values
        for (int i = 0; i < 2; i++) {
            gtk_entry_set_text(GTK_ENTRY(entries_line[i]), "");
        }
        gtk_label_set_text(GTK_LABEL(label_result_line), "");
        gtk_label_set_text(GTK_LABEL(label_m_line), "");
        gtk_label_set_text(GTK_LABEL(label_b_line), "");
        gtk_button_set_label(button, "Show line");
    } else {
        // Check if entries are valid
        for (int i = 0; i < 2; i++) {
            if (!GTK_IS_ENTRY(entries_line[i])) {
                //g_print("Invalid GtkEntry pointer at index %d\n", i);
                gtk_button_set_label(button, "Clear");
                return;
            }
        }

        const gchar *m_value = gtk_entry_get_text(GTK_ENTRY(entries_line[0]));
        const gchar *b_value = gtk_entry_get_text(GTK_ENTRY(entries_line[1]));

        // Debug statements to check the entered values
        //g_print("m_value: %s\n", m_value);
        //g_print("b_value: %s\n", b_value);

        // Convert the m_value to a double
        res_m_line = g_ascii_strtod(m_value, NULL);

        // Convert the b_value to a double
        res_b_line = g_ascii_strtod(b_value, NULL);

        // Debug statements to check the converted values
        //g_print("Converted res_m_line: %f\n", res_m_line);
        //g_print("Converted res_b_line: %f\n", res_b_line);

        // Display the result
        gchar *result_text = g_strdup_printf("y = %.2fx + %.2f", res_m_line, res_b_line);
        gtk_label_set_text(GTK_LABEL(label_result_line), result_text);
        g_free(result_text);

        // Display m and b values
        gchar *m_text = g_strdup_printf("m = %.2f", res_m_line);
        gtk_label_set_text(GTK_LABEL(label_m_line), m_text);
        g_free(m_text);

        gchar *b_text = g_strdup_printf("b = %.2f", res_b_line);
        gtk_label_set_text(GTK_LABEL(label_b_line), b_text);
        g_free(b_text);

        gtk_button_set_label(button, "Clear");

        draw_line = TRUE; // Set flag to true
        //width / 2g_print("Debug: draw_line set to TRUE\n"); // Debug statement

        if (GTK_IS_WIDGET(calcdrawing_area)) {
            gtk_widget_queue_draw(calcdrawing_area); // Force redraw of the drawing area
        } else {
            //g_print("Error: drawing_area is not a valid GtkWidget\n");
        }
    }
    is_empty = !is_empty;
}
///END LINEAR EQUATION

//SQUARE EQUATION
void on_eq_second_button_clicked(GtkWidget *button, GtkWidget *eq_second_box) {
    if (gtk_widget_get_visible(eq_second_box)) {
        gtk_widget_hide(eq_second_box);
    } else {
        gtk_widget_show(eq_second_box);
    }
}

void on_show_eq_second_button_clicked(GtkButton *button, gpointer user_data) {
    static gboolean is_empty = FALSE;
    gpointer *data = (gpointer *)user_data;
    GtkWidget **entries_eqsecond = (GtkWidget **)data[0];
    GtkWidget *label_result_eqsecond = GTK_WIDGET(data[1]);
    GtkWidget *label_x1_eqsecond = GTK_WIDGET(data[2]);
    GtkWidget *label_x2_eqsecond = GTK_WIDGET(data[3]); 
    GtkWidget *label_delta_eqsecond = GTK_WIDGET(data[4]); 

    if (is_empty) {
        // Reset all values
        for (int i = 0; i < 3; i++) {
            gtk_entry_set_text(GTK_ENTRY(entries_eqsecond[i]), "");
        }
        gtk_label_set_text(GTK_LABEL(label_result_eqsecond), "");
        gtk_label_set_text(GTK_LABEL(label_delta_eqsecond), "");
        gtk_label_set_text(GTK_LABEL(label_x1_eqsecond), "");
        gtk_label_set_text(GTK_LABEL(label_x2_eqsecond), "");
        gtk_button_set_label(button, "Show equation");
        //g_print("Button label set to 'Show equation'\n");

    } else {
        // Check if entries are valid
        for (int i = 0; i < 3; i++) {
            if (!GTK_IS_ENTRY(entries_eqsecond[i])) {
                //g_print("Invalid GtkEntry pointer at index %d\n", i);
                return;
            }
        }

        const gchar *a_value = gtk_entry_get_text(GTK_ENTRY(entries_eqsecond[0]));
        const gchar *b_value = gtk_entry_get_text(GTK_ENTRY(entries_eqsecond[1]));
        const gchar *c_value = gtk_entry_get_text(GTK_ENTRY(entries_eqsecond[2]));

        gchar *result_text = g_strdup_printf("%sx² + %sx + %s = 0", a_value, b_value, c_value);
        gtk_label_set_text(GTK_LABEL(label_result_eqsecond), result_text);
        g_free(result_text);

        double a = atof(a_value);
        double b = atof(b_value);
        double c = atof(c_value);
        a_eq2 = a;
        b_eq2 = b;
        c_eq2 = c;
        // Check if 'a' is zero to avoid division by zero
        if (a == 0) {
            gtk_label_set_text(GTK_LABEL(label_result_eqsecond), "Coefficient 'a' cannot be zero.");
            return;
        }

        double delta = (b * b) - (4 * a * c);
        gchar *delta_text = g_strdup_printf("Δ = %.2f", delta);
        gtk_label_set_text(GTK_LABEL(label_delta_eqsecond), delta_text);
        g_free(delta_text);
        delta_eq2 = delta;
        // Check if delta is negative to avoid calculating square root of a negative number
        if (delta < 0) {
            gtk_label_set_text(GTK_LABEL(label_result_eqsecond), "Delta is negative, no real roots.");
            gtk_label_set_text(GTK_LABEL(label_x1_eqsecond), "");
            gtk_label_set_text(GTK_LABEL(label_x2_eqsecond), "");
            return;
        }

        double x1_eqsecond = (-b + sqrt(delta)) / (2 * a);
        gchar *x1_text = g_strdup_printf("x₁ = %.2f", x1_eqsecond);
        gtk_label_set_text(GTK_LABEL(label_x1_eqsecond), x1_text);
        g_free(x1_text);
        x1_eq2 = x1_eqsecond;
        
        double x2_eqsecond = (-b - sqrt(delta)) / (2 * a);
        gchar *x2_text = g_strdup_printf("x₂ = %.2f", x2_eqsecond);
        gtk_label_set_text(GTK_LABEL(label_x2_eqsecond), x2_text);
        g_free(x2_text);
        x2_eq2 = x2_eqsecond;
        gtk_button_set_label(button, "Clear");
        //g_print("Button label set to 'Clear'\n");

        draw_eq2 = TRUE; // Set flag to true
        //g_print("Debug: draw_eq2 set to TRUE\n"); // Debug statement

        if (GTK_IS_WIDGET(calcdrawing_area)) {
            gtk_widget_queue_draw(calcdrawing_area); // Force redraw of the drawing area
        } else {
            //g_print("Error: drawing_area is not a valid GtkWidget\n");
        }
    }

    is_empty = !is_empty;
    //g_print("is_empty toggled, new value: %d\n", is_empty);
}

//END SQUARE EQUATION

//CUBIC EQUATION
void on_eq_cubic_button_clicked(GtkWidget *button, GtkWidget *eq_cubic_box) {
    if (gtk_widget_get_visible(eq_cubic_box)) {
        gtk_widget_hide(eq_cubic_box);
    } else {
        gtk_widget_show(eq_cubic_box);
    }
}

void on_show_eq_cubic_button_clicked(GtkButton *button, gpointer user_data) {
    static gboolean is_empty = FALSE;
    gpointer *data = (gpointer *)user_data;
    GtkWidget **entries_eqcubic = (GtkWidget **)data[0];
    GtkWidget *label_result_eqcubic = GTK_WIDGET(data[1]);
    GtkWidget *label_x_eqcubic= GTK_WIDGET(data[2]);


    // Check for null pointers
    if (!entries_eqcubic[0] || !entries_eqcubic[1] || !entries_eqcubic[2] || !entries_eqcubic[3]) {
        //g_print("Error: One or more entries_eqcubic pointers are null\n");
        return;
    }
    if (!label_result_eqcubic || !label_x_eqcubic) {
        //g_print("Error: One or more label pointers are null\n");
        return;
    }

    if (is_empty) {
        // Reset all values
        for (int i = 0; i < 4; i++) {
            gtk_entry_set_text(GTK_ENTRY(entries_eqcubic[i]), "");
        }
        gtk_label_set_text(GTK_LABEL(label_result_eqcubic), "");
        gtk_label_set_text(GTK_LABEL(label_x_eqcubic), "");
        gtk_button_set_label(button, "Show equation");
     } else {
        // Check if entries are valid
        for (int i = 0; i < 4; i++) {
            if (!GTK_IS_ENTRY(entries_eqcubic[i])) {
                //g_print("Invalid GtkEntry pointer at index %d\n", i);
                return;
            }
        }

        const gchar *a_value = gtk_entry_get_text(GTK_ENTRY(entries_eqcubic[0]));
        const gchar *b_value = gtk_entry_get_text(GTK_ENTRY(entries_eqcubic[1]));
        const gchar *c_value = gtk_entry_get_text(GTK_ENTRY(entries_eqcubic[2]));
        const gchar *d_value = gtk_entry_get_text(GTK_ENTRY(entries_eqcubic[3]));
        
        gchar *result_text = g_strdup_printf("%sx³ + %sx² + %sx + %s = 0", a_value, b_value, c_value, d_value);
        gtk_label_set_text(GTK_LABEL(label_result_eqcubic), result_text);
        g_free(result_text);

        double a = atof(a_value);
        double b = atof(b_value);
        double c = atof(c_value);
        double d = atof(d_value);
        a_eq3 = a;
        b_eq3 = b;
        c_eq3 = c;
        d_eq3 = d;

        // Check if 'a' is zero to avoid division by zero
        if (a == 0) {
            gtk_label_set_text(GTK_LABEL(label_result_eqcubic), "Coefficient 'a' cannot be zero.");
            return;
        }

        // Using Cardano's formula to find the real root
        double p = (3*a*c - b*b) / (3*a*a);
        double q = (2*b*b*b - 9*a*b*c + 27*a*a*d) / (27*a*a*a);
        double discriminant = (q*q) / 4 + (p*p*p) / 27;

        if (discriminant > 0) {
            double u = cbrt(-q/2 + sqrt(discriminant));
            double v = cbrt(-q/2 - sqrt(discriminant));
            double x1_eqcubic = u + v - b/(3*a);
            gchar *x_text = g_strdup_printf("x = %.2f", x1_eqcubic);
            gtk_label_set_text(GTK_LABEL(label_x_eqcubic), x_text);
            g_free(x_text);
            x_eq3 = x1_eqcubic;
        } else {
            gtk_label_set_text(GTK_LABEL(label_x_eqcubic), "No real roots.");
        }

        gtk_button_set_label(button, "Clear");
        draw_eq3 = TRUE; // Set flag to true

        if (GTK_IS_WIDGET(calcdrawing_area)) {
            gtk_widget_queue_draw(calcdrawing_area); // Force redraw of the drawing area
        } else {
            //g_print("Error: drawing_area is not a valid GtkWidget\n");
        }
    }

    is_empty = !is_empty;

}

//END CUBIC EQUATION

//ARCHIMEDEAN SPIRAL
void on_eq_archim_button_clicked(GtkWidget *button, GtkWidget *eq_archim_box) {
    if (gtk_widget_get_visible(eq_archim_box)) {
        gtk_widget_hide(eq_archim_box);
    } else {
        gtk_widget_show(eq_archim_box);
    }
}

void on_show_archim_button_clicked(GtkButton *button, gpointer user_data) {
    static gboolean is_empty = FALSE;
    gpointer *data = (gpointer *)user_data;
    GtkWidget **entries_eqarchim = (GtkWidget **)data[0];

    // Check for null pointers
    if (!entries_eqarchim[0] || !entries_eqarchim[1] || !entries_eqarchim[2]) {
        //g_print("Error: One or more entries_eqarchim pointers are null\n");
        return;
    }

    if (is_empty) {
        // Reset all values
        for (int i = 0; i < 3; i++) {
            gtk_entry_set_text(GTK_ENTRY(entries_eqarchim[i]), "");
        }
        gtk_button_set_label(button, "Show spiral");
    } else {
        // Check if entries are valid
        for (int i = 0; i < 3; i++) {
            if (!GTK_IS_ENTRY(entries_eqarchim[i])) {
                //g_print("Invalid GtkEntry pointer at index %d\n", i);
                return;
            }
        }

        const gchar *t_value = gtk_entry_get_text(GTK_ENTRY(entries_eqarchim[0]));
        const gchar *v_value = gtk_entry_get_text(GTK_ENTRY(entries_eqarchim[1]));
        const gchar *w_value = gtk_entry_get_text(GTK_ENTRY(entries_eqarchim[2]));

        //g_print("Text values - t: %s, v: %s, w: %s\n", t_value, v_value, w_value);

        double t = atof(t_value);
        double v = atof(v_value);
        double w = atof(w_value);
        
        // Store the values for later use
        t_eqarchim = t;
        v_eqarchim = v;
        w_eqarchim = w;

        gtk_button_set_label(button, "Clear");
                gtk_button_set_label(button, "Clear");
        draw_archim = TRUE; // Set flag to true

        if (GTK_IS_WIDGET(calcdrawing_area)) {
            gtk_widget_queue_draw(calcdrawing_area); // Force redraw of the drawing area
        } else {
            //g_print("Error: drawing_area is not a valid GtkWidget\n");
        }
    }

    is_empty = !is_empty;
}

//END ARCHIMEDEAN SPIRAL

//START EXPONENTIAL FUNCTION
void on_eq_expo_button_clicked(GtkWidget *button, GtkWidget *expo_box) {
    if (gtk_widget_get_visible(expo_box)) {
        gtk_widget_hide(expo_box);
    } else {
        gtk_widget_show(expo_box);
    }
}

void on_show_expo_button_clicked(GtkButton *button, gpointer user_data) {
    static gboolean is_empty = FALSE;
    gpointer *data = (gpointer *)user_data;
    GtkWidget **entries_expo = (GtkWidget **)data[0];

    // Check for null pointers
    if (!entries_expo[0] || !entries_expo[1] || !entries_expo[2] || !entries_expo[3] || !entries_expo[3] || !entries_expo[5] || !entries_expo[6]) {
        //g_print("Error: One or more entries_expo pointers are null\n");
        return;
    }

if (is_empty) {
    // Reset all values
    for (int i = 0; i < 7; i++) {
        gtk_entry_set_text(GTK_ENTRY(entries_expo[i]), "");
    }
    gtk_button_set_label(button, "Show function");
} else {
    // Check if entries are valid
    for (int i = 0; i < 7; i++) {
        if (!GTK_IS_ENTRY(entries_expo[i])) {
            //g_print("Invalid GtkEntry pointer at index %d\n", i);
            return;
        }
    }

    const gchar *x1_value = gtk_entry_get_text(GTK_ENTRY(entries_expo[0]));
    const gchar *x2_value = gtk_entry_get_text(GTK_ENTRY(entries_expo[1]));
    const gchar *x3_value = gtk_entry_get_text(GTK_ENTRY(entries_expo[2]));
    const gchar *x4_value = gtk_entry_get_text(GTK_ENTRY(entries_expo[3]));
    const gchar *x5_value = gtk_entry_get_text(GTK_ENTRY(entries_expo[4]));
    const gchar *x6_value = gtk_entry_get_text(GTK_ENTRY(entries_expo[5]));
    const gchar *b_value = gtk_entry_get_text(GTK_ENTRY(entries_expo[6]));       

    double x1 = atof(x1_value);
    double x2 = atof(x2_value);
    double x3 = atof(x3_value);
    double x4 = atof(x4_value);
    double x5 = atof(x5_value);
    double x6 = atof(x6_value);
    double b = atof(b_value);

    // Store the values for later use
    x1_expo = x1;
    x2_expo = x2;
    x3_expo = x3;
    x4_expo = x4;
    x5_expo = x5;
    x6_expo = x6;
    b_expo = b;

    if (b > 0) {
        draw_expo = TRUE; // Set flag to true only if b > 0
    } else {
        //g_print("Warning: b must be greater than 0 to set draw_expo to TRUE\n");
    }

    gtk_button_set_label(button, "Clear");

    if (GTK_IS_WIDGET(calcdrawing_area)) {
        gtk_widget_queue_draw(calcdrawing_area); // Force redraw of the drawing area
    } else {
        //g_print("Error: drawing_area is not a valid GtkWidget\n");
    }
}

    is_empty = !is_empty;
}
//END EXPONENTIAL FUNCTION

//START PROPORTIONS FUNCTION
void on_eq_proportions_button_clicked(GtkWidget *button, GtkWidget *proportions_box) {
    if (gtk_widget_get_visible(proportions_box)) {
        gtk_widget_hide(proportions_box);
    } else {
        gtk_widget_show(proportions_box);
    }
}

void on_show_proportions_button_clicked(GtkButton *button, gpointer user_data) {
    static gboolean is_empty = FALSE;
    gpointer *data = (gpointer *)user_data;
    GtkWidget **entries_proportions = (GtkWidget **)data[0];
    GtkWidget *label_result_proportions = GTK_WIDGET(data[1]);
    GtkWidget *label_x_proportions = GTK_WIDGET(data[2]);

    if (is_empty) {
        // Reset all values
        for (int i = 0; i < 4; i++) {
            gtk_entry_set_text(GTK_ENTRY(entries_proportions[i]), "");
        }
        gtk_label_set_text(GTK_LABEL(label_result_proportions), "");
        gtk_label_set_text(GTK_LABEL(label_x_proportions), "");
        gtk_button_set_label(button, "Show result");
        //g_print("Button label set to 'Show result'\n");
    } else {
        // Check if entries are valid
        for (int i = 0; i < 3; i++) {
            if (!GTK_IS_ENTRY(entries_proportions[i])) {
             //   g_print("Invalid GtkEntry pointer at index %d\n", i);
                return;
            }
        }

        const gchar *a_value = gtk_entry_get_text(GTK_ENTRY(entries_proportions[0]));
        const gchar *b_value = gtk_entry_get_text(GTK_ENTRY(entries_proportions[1]));
        const gchar *c_value = gtk_entry_get_text(GTK_ENTRY(entries_proportions[2]));
        const gchar *d_value = gtk_entry_get_text(GTK_ENTRY(entries_proportions[3]));  // Corrected the missing d_value

        gchar *result_text = g_strdup_printf("%s : %s = %s : %s", a_value, b_value, c_value, d_value);
        gtk_label_set_text(GTK_LABEL(label_result_proportions), result_text);
        g_free(result_text);

        double a = atof(a_value);
        double b = atof(b_value);
        double c = atof(c_value);
        double d = atof(d_value);  // Corrected the assignment for d
        double x = 0;

        // Check if 'a' is zero to avoid division by zero
        if (a == 0) {
            x = (b * c) / d;
            a = x;
        } 
        if (b == 0) {
            x = (d * a) / c;
            b = x;
        } 
        if (c == 0) {
            x = (d * a) / b;
            c = x;
        } 
        if (d == 0) {
            x = (c * b) / a;
            d = x;
        } 
        else {
        //    g_print("Error: 'a' value is zero, cannot calculate 'x'\n");
        }

        a_proportions = a;
        b_proportions = b;
        c_proportions = c;
        d_proportions = d;
        
         double highest_value = a;
        if (b > highest_value) highest_value = b;
        if (c > highest_value) highest_value = c;
        if (d > highest_value) highest_value = d;
   expression_result = highest_value;
   
        gchar *x_text = g_strdup_printf("x = %.2f", x);
        gtk_label_set_text(GTK_LABEL(label_x_proportions), x_text);
        g_free(x_text);

        gtk_button_set_label(button, "Clear");
        draw_proportions = TRUE; // Set flag to true
      //  g_print("Debug: draw_proportions set to TRUE\n"); // Debug statement

        if (GTK_IS_WIDGET(calcdrawing_area)) {
            gtk_widget_queue_draw(calcdrawing_area); // Force redraw of the drawing area
        } else {
    //        g_print("Error: calcdrawing_area is not a valid GtkWidget\n");
        }
    }

    is_empty = !is_empty;
  //  g_print("is_empty toggled, new value: %d\n", is_empty);
}

//END PROPORTION FUNCTION

//INEQ 2
void on_eq_ineq2_button_clicked(GtkWidget *button, GtkWidget *ineq2_box) {
    if (gtk_widget_get_visible(ineq2_box)) {
        gtk_widget_hide(ineq2_box);
    } else {
        gtk_widget_show(ineq2_box);
    }
}

void on_show_ineq2_button_clicked(GtkButton *button, gpointer user_data) {
    static gboolean is_empty = FALSE;
    gpointer *data_ineq2 = (gpointer *)user_data;
    GtkWidget **entries_ineq2 = (GtkWidget **)data_ineq2[0];
    GtkWidget *label_result_ineq2 = GTK_WIDGET(data_ineq2[1]);
    GtkWidget *label_x_ineq2 = GTK_WIDGET(data_ineq2[2]);
    GtkWidget *label_y1_ineq2 = GTK_WIDGET(data_ineq2[3]); 
    GtkWidget *label_x2_ineq2 = GTK_WIDGET(data_ineq2[4]);
    GtkWidget *label_y2_ineq2 = GTK_WIDGET(data_ineq2[5]);

    if (is_empty) {
        // Reset all values
        for (int i = 0; i < 3; i++) {
            gtk_entry_set_text(GTK_ENTRY(entries_ineq2[i]), "");
            }
        gtk_label_set_text(GTK_LABEL(label_result_ineq2), "");
        gtk_label_set_text(GTK_LABEL(label_x_ineq2), "");
        gtk_label_set_text(GTK_LABEL(label_y1_ineq2), "");
        gtk_label_set_text(GTK_LABEL(label_x2_ineq2), "");
        gtk_label_set_text(GTK_LABEL(label_y2_ineq2), "");
        gtk_button_set_label(button, "Show equation");
        //g_print("Button label set to 'Show equation'\n");
    } else {
        // Check if entries are valid
        for (int i = 0; i < 3; i++) {
            if (!GTK_IS_ENTRY(entries_ineq2[i])) {
                //g_print("Invalid GtkEntry pointer at index %d\n", i);
                return;
            }
        }

        const gchar *a_value = gtk_entry_get_text(GTK_ENTRY(entries_ineq2[0]));
        const gchar *b_value = gtk_entry_get_text(GTK_ENTRY(entries_ineq2[1]));
        const gchar *c_value = gtk_entry_get_text(GTK_ENTRY(entries_ineq2[2]));

        gchar *result_text = g_strdup_printf("%sx + %sy = %s", a_value, b_value, c_value);
        gtk_label_set_text(GTK_LABEL(label_result_ineq2), result_text);
        g_free(result_text);
   

        double a = atof(a_value);
        double b = atof(b_value);
        double c = atof(c_value);
        a_ineq2 = a;
        b_ineq2 = b;
        c_ineq2 = c;

        // Check if 'a' is zero to avoid division by zero
        if (a == 0) {
            gtk_label_set_text(GTK_LABEL(label_result_ineq2), "Coefficient 'a' cannot be zero.");
            return;
        }

    x1_ineq2 = 0;
    y1_ineq2 = c_ineq2 / b_ineq2;
    gchar *x1_text = g_strdup_printf("x1 = %.3f", x1_ineq2);
    gchar *y1_text = g_strdup_printf("y1 = %.3f", y1_ineq2);
    gtk_label_set_text(GTK_LABEL(label_x_ineq2), x1_text);
    gtk_label_set_text(GTK_LABEL(label_y1_ineq2), y1_text);
    g_free(x1_text);
    g_free(y1_text);

 
    x2_ineq2 = c_ineq2;
    double ax = (a_ineq2 * x2_ineq2);
    
    if (b_ineq2 > 0) {
    c_ineq2 = c_ineq2 - ax;
    y2_ineq2 = c_ineq2 / b;   
    gchar *x2_text = g_strdup_printf("x2 = %.3f", x2_ineq2);
    gchar *y2_text = g_strdup_printf("y2 = %.3f", y2_ineq2);
    gtk_label_set_text(GTK_LABEL(label_x2_ineq2), x2_text);
    gtk_label_set_text(GTK_LABEL(label_y2_ineq2), y2_text);
    g_free(x2_text);
    g_free(y2_text);
    }
    
    if (b_ineq2 < 0) {
    if (a > 0) {
    c_ineq2 = c_ineq2 - ax;
    }
    if (a < 0)  {
    c_ineq2 = c_ineq2 - ax;
    }
    //g_print("ax = %.3f\n", ax);
    //g_print("c_ineq2 = %.3f\n", c_ineq2);
    y2_ineq2 = (c_ineq2 / b);   
    gchar *x2_text = g_strdup_printf("x2 = %.3f", x2_ineq2);
    gchar *y2_text = g_strdup_printf("y2 = %.3f", y2_ineq2);
    gtk_label_set_text(GTK_LABEL(label_x2_ineq2), x2_text);
    gtk_label_set_text(GTK_LABEL(label_y2_ineq2), y2_text);
    g_free(x2_text);
    g_free(y2_text);
    }
    


        gtk_button_set_label(button, "Clear");
 
        draw_ineq2 = TRUE; // Set flag to true
      //  g_print("Debug: draw_proportions set to TRUE\n"); // Debug statement

        if (GTK_IS_WIDGET(calcdrawing_area)) {
            gtk_widget_queue_draw(calcdrawing_area); // Force redraw of the drawing area
        } else {
    //        g_print("Error: calcdrawing_area is not a valid GtkWidget\n");
        }
    }

    is_empty = !is_empty;
  //  g_print("is_empty toggled, new value: %d\n", is_empty);
}
/*end calc*/

static WebKitWebView *web_view;
static WebKitWebView *hidden_web_view;
static void download_started_callback(WebKitDownload *download, gpointer user_data);
static void load_finished(WebKitDownload *download, gpointer user_data);
static void progress_updated(WebKitDownload *download, GParamSpec *pspec, gpointer user_data);

void enable_experimental_features(WebKitWebView* web_view) {
    if (!web_view) return;

    WebKitSettings* settings = webkit_web_view_get_settings(web_view);
    webkit_settings_set_enable_media_stream(settings, TRUE);
    webkit_settings_set_enable_webgl(settings, TRUE);
    webkit_settings_set_enable_media_capabilities(settings, TRUE);
    webkit_settings_set_enable_site_specific_quirks(settings, TRUE);
    webkit_settings_set_enable_developer_extras(settings, TRUE);    
}



//global declaration
GtkWidget *entry;
GtkWidget *download_event_box;
GtkWidget *downloadbar;
void download_started_callback(WebKitDownload *download, gpointer user_data);



WebKitWebView* setup_webview_storage() {
    // Define storage paths
    const gchar *base_directory = "/usr/share/debris/main_file/";
    const gchar *cache_directory = g_build_filename(base_directory, "caches.txt", NULL);
    const gchar *cookie_directory = g_build_filename(base_directory, "cookies.txt", NULL);

    // Create a Website Data Manager with cache and cookie storage paths
    WebKitWebsiteDataManager *data_manager = webkit_website_data_manager_new(
        nullptr,  // Local storage path
        cache_directory,  // Disk cache path
        nullptr,  // IndexedDB path
        nullptr,  // WebSQL path
        nullptr,  // Cache storage path
        cookie_directory  // Cookie storage path
    );

// Create a Web Context and set the data manager
WebKitWebContext *web_context = webkit_web_context_new_with_website_data_manager(data_manager);

// Connect the download signal AFTER creating `web_context`
g_signal_connect(web_context, "download-started", G_CALLBACK(download_started_callback), NULL);

// Create WebView with the Web Context
WebKitWebView *web_view = WEBKIT_WEB_VIEW(webkit_web_view_new_with_context(web_context));


    // Get the Cookie Manager and enable persistent storage
    WebKitCookieManager *cookie_manager = webkit_web_context_get_cookie_manager(web_context);
    webkit_cookie_manager_set_persistent_storage(cookie_manager, cookie_directory, WEBKIT_COOKIE_PERSISTENT_STORAGE_SQLITE);

    // Set cookie accept policy
    webkit_cookie_manager_set_accept_policy(cookie_manager, WEBKIT_COOKIE_POLICY_ACCEPT_NO_THIRD_PARTY);


    WebKitSettings *settings = webkit_web_view_get_settings(web_view);
    g_object_set(G_OBJECT(settings),
        "enable-webgl", FALSE,  // Disable WebGL (potential tracking vector)
        "enable-java", FALSE,  // Disable Java applets
        "enable-javascript", TRUE,  // Keep JavaScript enabled (can be disabled if needed)
        NULL);

    return web_view;
}



// Function to load CSS file
static void load_css(const char *css_file) {
    GtkCssProvider *provider = gtk_css_provider_new();
    gboolean success = gtk_css_provider_load_from_path(provider, css_file, NULL);
    if (success) {
        gtk_style_context_add_provider_for_screen(gdk_screen_get_default(),
                                                  GTK_STYLE_PROVIDER(provider),
                                                  GTK_STYLE_PROVIDER_PRIORITY_USER);
        //g_print("CSS file '%s' loaded successfully.\n", css_file);
    } else {
        //g_print("Failed to load CSS file '%s'.\n", css_file);
    }
    g_object_unref(provider);
}



//SAVE PREFERENCES THEME

void save_preference(const char *key, const char *value) {
    GKeyFile *key_file = g_key_file_new();
    char *config_path = g_build_filename(g_get_user_config_dir(), "debris_config.ini", NULL);

    // Ensure the directory exists (important for user-specific config)
    g_mkdir_with_parents(g_get_user_config_dir(), 0755);

    // Check if file exists, and initialize it if missing
    struct stat buffer;
    if (stat(config_path, &buffer) != 0) {
        printf("Config file not found. Creating new file: %s\n", config_path);

        // Create an empty file before writing preferences
        FILE *fp = fopen(config_path, "w");
        if (fp) {
            fclose(fp);
        } else {
            fprintf(stderr, "Error: Unable to create config file.\n");
            g_free(config_path);
            g_key_file_free(key_file);
            return;
        }
    }

    // Load existing config file or continue with a fresh instance
    GError *error = NULL;
    g_key_file_load_from_file(key_file, config_path, G_KEY_FILE_NONE, NULL);

    // Update the preference
    g_key_file_set_string(key_file, "Preferences", key, value);

    // Save the updated config file
    if (!g_key_file_save_to_file(key_file, config_path, &error)) {
        fprintf(stderr, "Error saving preference [%s]: %s\n", key, error->message);
        g_error_free(error);
    } else {
        printf("Preference saved [%s]: %s\n", key, value);
    }

    g_free(config_path);
    g_key_file_free(key_file);
}

// Saves selected theme
void save_theme_choice(const char *theme) {
    save_preference("theme", theme);
}

//THEMES
static void on_themebutton_light_clicked(GtkButton *themebutton_light, gpointer user_data) {
    load_css("/usr/share/debris/css/light_theme.css"); // Load css2.css when the button is clicked
    save_theme_choice("/usr/share/debris/css/light_theme.css"); 
}

static void on_themebutton_dark_clicked(GtkButton *themebutton_dark, gpointer user_data) {
    load_css("/usr/share/debris/css/dark_theme.css"); // Load css2.css when the button is clicked
    save_theme_choice("/usr/share/debris/css/dark_theme.css"); 
}

static void on_themebutton_darkcolor_clicked(GtkButton *themebutton_darkcolor, gpointer user_data) {
    load_css("/usr/share/debris/css/darkcolor_theme.css"); // Load css2.css when the button is clicked
    save_theme_choice("/usr/share/debris/css/darkcolor_theme.css"); 
}

static void on_themebutton_lightcolor_clicked(GtkButton *themebutton_lightcolor, gpointer user_data) {
    load_css("/usr/share/debris/css/lightcolor_theme.css"); // Load css2.css when the button is clicked
    save_theme_choice("/usr/share/debris/css/lightcolor_theme.css");
}

static void on_themebutton_matrix_clicked(GtkButton *themebutton_matrix, gpointer user_data) {
    load_css("/usr/share/debris/css/matrix_theme.css"); // Load css2.css when the button is clicked
    save_theme_choice("/usr/share/debris/css/matrix_theme.css");
}

static void on_themebutton_accessibility1_clicked(GtkButton *themebutton_accessibility1, gpointer user_data) {
    load_css("/usr/share/debris/css/accessibility1.css"); // Load css2.css when the button is clicked
    save_theme_choice("/usr/share/debris/css/accessibility1.css");
}

static void on_themebutton_accessibility_dark_clicked(GtkButton *themebutton_accessibility_dark, gpointer user_data) {
    load_css("/usr/share/debris/css/accessibility_dark.css"); // Load css2.css when the button is clicked
    save_theme_choice("/usr/share/debris/css/accessibility_dark.css");
}

static void on_themebutton_accessibility_dark_cyan_clicked(GtkButton *themebutton_accessibility_dark_cyan, gpointer user_data) {
    load_css("/usr/share/debris/css/accessibility_dark_cyan.css"); // Load css2.css when the button is clicked
    save_theme_choice("/usr/share/debris/css/accessibility_dark_cyan.css");
}

//load saved themes
void load_saved_theme() {
    GKeyFile *key_file = g_key_file_new();
    char *config_path = g_build_filename(g_get_user_config_dir(), "debris_config.ini", NULL);

    GError *error = NULL;
    if (!g_key_file_load_from_file(key_file, config_path, G_KEY_FILE_NONE, &error)) {
        fprintf(stderr, "Error loading config file: %s\n", error->message);
        g_error_free(error);

        // Ensure the config file exists
        g_key_file_set_string(key_file, "Preferences", "theme", "/usr/share/debris/css/dark_theme.css");
        g_key_file_save_to_file(key_file, config_path, NULL);

        load_css("/usr/share/debris/css/dark_theme.css"); // Default to dark theme
        g_free(config_path);
        g_key_file_free(key_file);
        return;
    }

    printf("Successfully loaded config file: %s\n", config_path);

    char *theme = g_key_file_get_string(key_file, "Preferences", "theme", NULL);
    if (theme && g_file_test(theme, G_FILE_TEST_EXISTS)) {
        load_css(theme); // Load the saved theme
        printf("Loaded theme: %s\n", theme);
        g_free(theme);
    } else {
        load_css("/usr/share/debris/css/dark_theme.css"); // Default to dark theme
        printf("Loaded default theme: dark_theme.css\n");
    }

    g_free(config_path);
    g_key_file_free(key_file);
}

//END THEMES

void on_uri_changed(WebKitWebView *web_view, GParamSpec *pspec, gpointer data);
void on_url_entry_activate(GtkEntry *entry, gpointer user_data);
void on_switch_page(GtkNotebook *notebook, GtkWidget *page, guint page_num, gpointer data);

static gboolean on_buttonmove_press_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
    if (event->button == GDK_BUTTON_PRIMARY) {
        gtk_window_begin_move_drag(GTK_WINDOW(widget), event->button, event->x_root, event->y_root, event->time);
    } else if (event->button == GDK_BUTTON_SECONDARY) {
        GdkWindowEdge edge;
        
        if (event->x < 10) {
            edge = GDK_WINDOW_EDGE_WEST;
        } else if (event->x > gtk_widget_get_allocated_width(widget) - 10) {
            edge = GDK_WINDOW_EDGE_EAST;
        } else if (event->y < 10) {
            edge = GDK_WINDOW_EDGE_NORTH;
        } else if (event->y > gtk_widget_get_allocated_height(widget) - 10) {
            edge = GDK_WINDOW_EDGE_SOUTH;
        } else {
            edge = GDK_WINDOW_EDGE_SOUTH_EAST;
        }
        
        gtk_window_begin_resize_drag(GTK_WINDOW(widget), edge, event->button, event->x_root, event->y_root, event->time);
    }
    return TRUE;
}

//BACK
void on_back_button_clicked(GtkButton *button, gpointer user_data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(user_data);
    int page_num = gtk_notebook_get_current_page(notebook);
    GtkWidget *current_page = gtk_notebook_get_nth_page(notebook, page_num);
    WebKitWebView *web_view = WEBKIT_WEB_VIEW(current_page);
    if (webkit_web_view_can_go_back(web_view)) {
        webkit_web_view_go_back(web_view);
    } else {
       // g_print("Cannot go back, no previous page available.\n");
    }
}

//
// FORWARD
void on_forward_button_clicked(GtkButton *button, gpointer user_data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(user_data);
    int page_num = gtk_notebook_get_current_page(notebook);
    GtkWidget *current_page = gtk_notebook_get_nth_page(notebook, page_num);
    WebKitWebView *web_view = WEBKIT_WEB_VIEW(current_page);
    if (webkit_web_view_can_go_forward(web_view)) {
        webkit_web_view_go_forward(web_view);
    } else {
       // g_print("Cannot go forward, no next page available.\n");
    }
}


// Callback function to update the tab label when the title changes
void on_title_changed(WebKitWebView *web_view, GParamSpec *pspec, gpointer user_data);

// Callback function to handle the creation of a new page
void on_new_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    // Create a new web view
    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://duckduckgo.com/");
    // Create a label for the tab
    GtkWidget *tab_label = gtk_label_new("New Tab");
    // Add the new web view to the notebook with the tab label
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    // Connect the title-changed signal to update the tab label
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), data);

 
    // Show all widgets in the new tab
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

///////
static void load_failed(WebKitWebView *web_view, WebKitLoadEvent load_event, const gchar *failing_uri, GError *error, gpointer user_data) {
    // Set the entry text to "error"
    gchar *current_dir = g_get_current_dir();
    gchar *error_file_path = g_strconcat("file://", current_dir, "/usr/share/debris/main_file/error.html", NULL);
    webkit_web_view_load_uri(web_view, error_file_path);
    g_free(current_dir);
    g_free(error_file_path);
}

// update the tab label when the title changes
void on_title_changed(WebKitWebView *web_view, GParamSpec *pspec, gpointer user_data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(gtk_widget_get_parent(GTK_WIDGET(web_view)));
    int page_num = GPOINTER_TO_INT(user_data);

    const gchar *title = webkit_web_view_get_title(web_view);
    if (title) {
        // Validate if the title is a valid UTF-8 string
        if (!g_utf8_validate(title, -1, NULL)) {
            g_warning("Invalid UTF-8 string detected. Using a fallback title.");
            title = "Invalid Title"; // Fallback title
        }

        gchar *truncated_title = NULL;

        // Check if the address starts with "https://www.google.com/search?q="
        const gchar *google_prefix = "https://www.google.com/search?q=";
        if (g_str_has_prefix(title, google_prefix)) {
            // Extract the part between "=" and "&"
            const gchar *query_start = strstr(title, "=");
            if (query_start) {
                query_start++; // Move past '='
                const gchar *query_end = strstr(query_start, "&");
                if (query_end) {
                    // Extract content between '=' and '&'
                    truncated_title = g_strndup(query_start, query_end - query_start);
                } else {
                    // No '&', take the rest of the string after '='
                    truncated_title = g_strdup(query_start);
                }
            }
        } else {
            // Fallback to truncating the title if not a Google search URL
            truncated_title = g_utf8_strncpy((gchar *)g_malloc0(11), title, 10);
        }

        if (truncated_title) {
            // Set the extracted or truncated title as the tab label
            GtkWidget *tab_label = gtk_notebook_get_tab_label(notebook, GTK_WIDGET(web_view));
            gtk_label_set_text(GTK_LABEL(tab_label), truncated_title);
            g_free(truncated_title);
        }
    }
}


void update_entry_with_current_url(WebKitWebView *web_view, GtkEntry *entry) {
    const gchar *uri = webkit_web_view_get_uri(web_view);

    // Print the URL to the console
    g_print("Current URL: %s\n", uri);

    if (GTK_IS_ENTRY(entry)) {
        gtk_entry_set_text(entry, uri);
        gtk_widget_queue_draw(GTK_WIDGET(entry)); // Force UI update
    } else {
        g_warning("Provided widget is not a GtkEntry");
    }
}

//fake switch to update entry correctly
void fake_tab_switch(GtkNotebook *notebook, GtkWidget *web_view) {
    // Get the current page number of the web view
    gint current_page = gtk_notebook_page_num(notebook, web_view);

    // Check if the web view is valid and part of the notebook
    if (current_page != -1) {
        // Temporarily set a different tab as active
        gint next_page = (current_page + 1) % gtk_notebook_get_n_pages(notebook);
        gtk_notebook_set_current_page(notebook, next_page);

        // Switch back to the original tab
        gtk_notebook_set_current_page(notebook, current_page);
    } else {
        g_warning("Web view not found in notebook");
    }
}

//HISTORY & ON URI_CHANGE
// Function to get the correct user-specific history path
std::string get_history_file_path() {
    const gchar *config_dir = g_get_user_config_dir(); // Detect user-specific config directory
    std::string debris_dir = std::string(config_dir) + "/debris";

    // Ensure the ~/.config/debris/ directory exists
    g_mkdir_with_parents(debris_dir.c_str(), 0755);

    return debris_dir + "/web_history.txt";
}

// Global history list
std::vector<std::pair<std::string, std::string>> history;

// Function declarations
std::string trim_url_history(const std::string& url);
std::string current_date_time();
void save_history();
void set_full_history_url(GtkWidget *button, const gchar *url);
const gchar *get_full_history_url(GtkWidget *button);

// Function to set notebook data for history
void set_history_notebook_data(GtkWidget *historybar, GtkNotebook *notebook) {
    g_object_set_data(G_OBJECT(historybar), "notebook", notebook);
    g_object_set_data(G_OBJECT(notebook), "historybar", historybar);

    // Initialize history on setup
    history.clear();
}



//Global function
void add_history_to_bar(GtkEntry *entry, GtkNotebook *notebook) ;

// Updated on_uri_changed function to save history
void on_uri_changed(WebKitWebView *web_view, GParamSpec *pspec, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(gtk_widget_get_parent(GTK_WIDGET(web_view)));
    GtkWidget *historybar = GTK_WIDGET(g_object_get_data(G_OBJECT(notebook), "historybar"));

    // Perform a fake tab switch
    fake_tab_switch(notebook, GTK_WIDGET(web_view));

    // Retrieve the current URI
    const gchar *uri = webkit_web_view_get_uri(web_view);
    if (uri) {
        //g_print("Current URI: %s\n", uri);

        // Trim the URI and save it with a timestamp
        std::string trimmed_uri = trim_url_history(uri);
        std::string timestamp = current_date_time();
        history.emplace_back(trimmed_uri, timestamp);
        save_history();

        // Add the new history item to the historybar dynamically
        if (historybar) {
            GtkEntry *dummy_entry = GTK_ENTRY(gtk_entry_new());
            gtk_entry_set_text(dummy_entry, uri);
            add_history_to_bar(dummy_entry, notebook);
            gtk_widget_destroy(GTK_WIDGET(dummy_entry)); // Clean up dummy entry
        }
    } else {
        g_print("No URI available.\n");
    }
}


// Trim the URL function (unchanged)
std::string trim_url_history(const std::string& url) {
    if (!url.empty() && url.back() == '/') {
        return url.substr(0, url.size() - 1); // Remove the trailing slash
    }
    return url;
}


// Save history function with duplicate URL check
void save_history() {
    std::ofstream file(get_history_file_path(), std::ios::out); // Explicitly open in output mode
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing: " << get_history_file_path() << std::endl;

        return;
    }

    if (history.empty()) {
        std::cout << "Debug: History is empty, nothing to save." << std::endl;
    } else {
        std::set<std::string> unique_urls; // Create a set to store unique URLs
        for (const auto &entry : history) {
            if (unique_urls.find(entry.first) == unique_urls.end()) { // Check if URL is already in the set
                file << entry.second << " | " << entry.first << "\n"; // Save timestamp and URL
                unique_urls.insert(entry.first); // Add URL to the set
            }
        }
    }

    file.close();
    if (file.fail()) {
        std::cerr << "Error occurred while closing the file." << std::endl;
    }
}

// Function to get the current date and time (unchanged)
std::string current_date_time() {
    time_t t = time(nullptr);
    tm *tm = localtime(&t);
    char buffer[80];
    if (strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", tm) == 0) {
        std::cerr << "Error formatting time.\n";
        return "Unknown Time";
    }
    return std::string(buffer);
}


// Function to remove a line containing a specific URL from web_history.txt
void remove_history_from_file(const gchar *url) {
    std::string history_file = get_history_file_path();
    std::string temp_file = history_file + ".temp";

    FILE *file = fopen(history_file.c_str(), "r");
    FILE *temp = fopen(temp_file.c_str(), "w");

    if (!file || !temp) {
        std::cerr << "Error opening history files." << std::endl;
        return;
    }

    gchar line[256];
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, url) == NULL) {
            fputs(line, temp);
        }
    }

    fclose(file);
    fclose(temp);
    remove(history_file.c_str());
    rename(temp_file.c_str(), history_file.c_str());

    std::cout << "History entry removed successfully." << std::endl;
}


void set_full_history_url(GtkWidget *button, const gchar *url) {
    g_object_set_data_full(G_OBJECT(button), "full-url", g_strdup(url), g_free);
}

const gchar *get_full_history_url(GtkWidget *button) {
    return (const gchar *)g_object_get_data(G_OBJECT(button), "full-url");
}


// Modified callback function to handle the creation of a new page with the history URL or deletion
void on_history_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    if (!GTK_IS_NOTEBOOK(notebook)) {
        std::cerr << "Error: Invalid GtkNotebook pointer." << std::endl;
        return;
    }

    const gchar *url = get_full_history_url(button);

// Create a new dialog
GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(gtk_widget_get_toplevel(button)),
                                           GTK_DIALOG_MODAL,
                                           GTK_MESSAGE_OTHER,
                                           GTK_BUTTONS_NONE,
                                           NULL);
gtk_widget_set_name(dialog, "dialog");

// Create a header bar
GtkWidget *header_bar = gtk_header_bar_new();
gtk_header_bar_set_show_close_button(GTK_HEADER_BAR(header_bar), TRUE);
gtk_header_bar_set_title(GTK_HEADER_BAR(header_bar), "Preview");
gtk_widget_set_name(header_bar, "dialog_header_bar");

// Add the header bar to the dialog's custom title
gtk_window_set_titlebar(GTK_WINDOW(dialog), header_bar);

// Create a container to hold the preview web view
GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

// Create a web view for the preview
WebKitWebView *preview_web_view = setup_webview_storage();
gtk_widget_set_size_request(GTK_WIDGET(preview_web_view), 300, 400);
webkit_web_view_load_uri(preview_web_view, url);

// Add the preview web view to the content area
gtk_box_pack_start(GTK_BOX(content_area), GTK_WIDGET(preview_web_view), TRUE, TRUE, 0);

// Add buttons to the dialog
GtkWidget *open_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Open", GTK_RESPONSE_ACCEPT);
gtk_widget_set_name(open_button, "no_button");
GtkWidget *delete_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Delete", GTK_RESPONSE_REJECT);
gtk_widget_set_name(delete_button, "yes_button");

// Show all the widgets
gtk_widget_show_all(dialog);

gint response = gtk_dialog_run(GTK_DIALOG(dialog));
gtk_widget_destroy(dialog);

// Handle dialog response
if (response == GTK_RESPONSE_ACCEPT) {
    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, url);
    enable_experimental_features(web_view);

    GtkWidget *tab_label = gtk_label_new(url);
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);

    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);

    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
} else if (response == GTK_RESPONSE_REJECT) {
    // Remove the history entry from the file
    gtk_widget_destroy(button);
    remove_history_from_file(get_history_file_path().c_str());

}

}



std::string wrap_text(const std::string &input, size_t width) {
    std::string output;
    size_t current_pos = 0;

    while (current_pos < input.length()) {
        size_t next_pos = current_pos + width;
        if (next_pos >= input.length()) {
            output += input.substr(current_pos);
            break;
        } else {
            size_t wrap_pos = next_pos;
            while (wrap_pos > current_pos && input[wrap_pos] != ' ') {
                wrap_pos--;
            }
            if (wrap_pos == current_pos) {
                wrap_pos = next_pos;
            }
            output += input.substr(current_pos, wrap_pos - current_pos) + '\n';
            current_pos = wrap_pos + 1;
        }
    }
    return output;
}

void load_history(GtkWidget *historybar) {
    std::ifstream file(get_history_file_path());
    if (!file) {
        std::cerr << "Error opening file for writing: " << get_history_file_path() << std::endl;

        return;
    }
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue; // Skip empty lines
        }

        // Find the separator " | "
        std::size_t separator_pos = line.find(" | ");
        if (separator_pos == std::string::npos) {
            std::cerr << "Error: Invalid format in history file: " << line << std::endl;
            continue;
        }

        std::string timestamp = line.substr(0, separator_pos);
        std::string url = line.substr(separator_pos + 3); // Skip over " | "

        if (url.empty() || timestamp.empty()) {
            std::cerr << "Error: Empty timestamp or URL in history file: " << line << std::endl;
            continue;
        }

        history.emplace_back(url, timestamp);

        // Get only the first 50 characters of the URL
        std::string shortened_url = url.substr(0, 50);

        // Create history label with shortened URL
        std::string history_label = "SITE: " + shortened_url + "\nDATE: " + timestamp;
        GtkWidget *history_button = gtk_button_new();
        GtkLabel *label = GTK_LABEL(gtk_label_new(history_label.c_str()));
        gtk_label_set_line_wrap(label, TRUE); // Enable line wrap for the label
        gtk_label_set_xalign(GTK_LABEL(label), 0.0); // Align text to the far left

        //gtk_widget_set_valign(GTK_WIDGET(label), GTK_ALIGN_START);
        gtk_container_add(GTK_CONTAINER(history_button), GTK_WIDGET(label));

        gtk_widget_set_name(history_button, "history_button");
        gtk_box_pack_start(GTK_BOX(historybar), history_button, FALSE, FALSE, 0); // Set expand to TRUE
        gtk_widget_set_halign(history_button, GTK_ALIGN_FILL); // Set horizontal alignment to FILL
        gtk_widget_show_all(history_button);

        set_full_history_url(history_button, url.c_str());

        GtkNotebook *notebook = GTK_NOTEBOOK(g_object_get_data(G_OBJECT(historybar), "notebook"));
        if (!GTK_IS_NOTEBOOK(notebook)) {
            std::cerr << "Error: Invalid GtkNotebook pointer in historybar." << std::endl;
            continue;
        }

        g_signal_connect(history_button, "clicked", G_CALLBACK(on_history_clicked), notebook);
    }
}



//remove double to solve google issue
// Preprocess URL to remove unnecessary query parameters
std::string preprocess_url(const std::string &url) {
    // Only target "google.com" URLs
    if (url.find("google.com") != std::string::npos) {
        size_t question_mark_pos = url.find('?');
        if (question_mark_pos != std::string::npos) {
            std::string base_url = url.substr(0, question_mark_pos); // Base URL before "?"
            std::string query_params = url.substr(question_mark_pos + 1); // Query parameters after "?"

            // Keep essential parameters and remove session-specific ones using regex
            std::regex tracking_params_regex("&(sca_esv|iflsig|sei)=[^&]*");
            query_params = std::regex_replace(query_params, tracking_params_regex, "");

            // Reconstruct the URL with cleaned query parameters
            return base_url + "?" + query_params;
        }
    }
    return url; // Return original URL if it's not from "google.com"
}

void remove_duplicates_from_history(GtkWidget *historybar) {
    std::string history_file = get_history_file_path(); // Get per-user file path

    std::ifstream file(history_file, std::ios::in); // Open the file for reading
    if (!file.is_open()) {
        std::cerr << "Error opening file for reading: " << history_file << std::endl;
        return;
    }

    // Use a set to store only unique preprocessed URLs
    std::set<std::string> unique_urls;
    std::vector<std::string> cleaned_lines; // Store cleaned lines for writing
    std::string line;

    // Read each line from the file
    while (std::getline(file, line)) {
        size_t separator_pos = line.find(" | ");
        if (separator_pos != std::string::npos) {
            std::string url = line.substr(separator_pos + 3); // Extract URL
            std::string preprocessed_url = preprocess_url(url); // Preprocess URL
            if (unique_urls.find(preprocessed_url) == unique_urls.end()) { // If preprocessed URL is unique
                unique_urls.insert(preprocessed_url); // Add preprocessed URL to the set
                cleaned_lines.push_back(line); // Keep original line (timestamp + URL) for writing back
            }
        } else {
            std::cerr << "Error: Line format invalid, skipping: " << line << std::endl;
        }
    }

    file.clear(); // Clear any error flags before closing the file
    file.close();
    if (file.fail()) {
        std::cerr << "Error occurred while closing the file after reading." << std::endl;
    }

    // Reopen the file for writing (overwriting the file)
    std::ofstream outfile(history_file, std::ios::out | std::ios::trunc);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file for writing: " << history_file << std::endl;
        return;
    }

    // Write back only cleaned lines
    for (const auto &cleaned_line : cleaned_lines) {
        outfile << cleaned_line << "\n";
    }

    outfile.close();
    if (outfile.fail()) {
        std::cerr << "Error occurred while closing the file after writing." << std::endl;
    }

    // Reload the history bar
    if (historybar) {
        load_history(historybar);
    } else {
        std::cerr << "Error: historybar is NULL, cannot reload history." << std::endl;
    }
}


void add_history_to_bar(GtkEntry *entry, GtkNotebook *notebook) {
    const gchar *url = gtk_entry_get_text(entry);
    if (url && *url) { // Ensure the URL is not NULL or empty
        std::string trimmed_uri = trim_url_history(url);
        std::string timestamp = current_date_time();

        // Wrap URL text after 50 characters
        std::string wrapped_uri = wrap_text(trimmed_uri, 50);

        std::string history_label = "SITE: " + wrapped_uri + "\nDATE: " + timestamp;

        // Add to the global history vector and save
        history.emplace_back(trimmed_uri, timestamp);
        save_history();

        // Retrieve the historybar from the notebook data
        GtkWidget *historybar = GTK_WIDGET(g_object_get_data(G_OBJECT(notebook), "historybar"));
        if (historybar) {
            // Remove duplicate entries from the history file
            remove_duplicates_from_history(historybar); // Now uses per-user history file

            // Clear current historybar widgets
            GList *children = gtk_container_get_children(GTK_CONTAINER(historybar));
            for (GList *child = children; child != NULL; child = child->next) {
                gtk_widget_destroy(GTK_WIDGET(child->data));
            }
            g_list_free(children);

            // Reload history from file
            load_history(historybar);

            // Force the historybar to update immediately
            gtk_widget_queue_draw(historybar);
        } else {
            std::cerr << "Error: History bar not found." << std::endl;
        }
    } else {
        std::cerr << "Error: No valid URL found in the GtkEntry." << std::endl;
    }
}


static void on_delall_history_button_clicked(GtkWidget *widget, gpointer data) {
    // Get the top-level window of the widget
    GtkWidget *toplevel = gtk_widget_get_toplevel(widget);

    if (GTK_IS_WIDGET(toplevel) && GTK_IS_WINDOW(toplevel)) {
        // Create a confirmation dialog
        GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(toplevel),
                                                   GTK_DIALOG_MODAL,
                                                   GTK_MESSAGE_QUESTION,
                                                   GTK_BUTTONS_NONE,
                                                   "Are you sure?");
        gtk_widget_set_name(dialog, "dialog");

        // Custom header bar for the dialog
        GtkWidget *header_bar = gtk_header_bar_new();
        gtk_header_bar_set_show_close_button(GTK_HEADER_BAR(header_bar), TRUE);
        gtk_header_bar_set_title(GTK_HEADER_BAR(header_bar), "Delete history");
        gtk_widget_set_name(header_bar, "dialog_header_bar");
        gtk_window_set_titlebar(GTK_WINDOW(dialog), header_bar);

        // Add dialog buttons
        GtkWidget *yes_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Yes", GTK_RESPONSE_ACCEPT);
        gtk_widget_set_name(yes_button, "yes_button");
        GtkWidget *no_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "No", GTK_RESPONSE_REJECT);
        gtk_widget_set_name(no_button, "no_button");

        gtk_widget_show_all(dialog);
        gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_CANCEL);

        // Handle the dialog response
        gint response = gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);

        if (response == GTK_RESPONSE_ACCEPT) {
            // Clear the user's history file
            FILE *file = fopen(get_history_file_path().c_str(), "w");
            if (file != NULL) {
                fclose(file);
            }

            // Reload the history UI dynamically
            GtkNotebook *notebook = GTK_NOTEBOOK(data);
            GtkWidget *historybar = GTK_WIDGET(g_object_get_data(G_OBJECT(notebook), "historybar"));

            if (historybar) {
                // Clear current historybar widgets
                GList *children = gtk_container_get_children(GTK_CONTAINER(historybar));
                for (GList *child = children; child != NULL; child = child->next) {
                    gtk_widget_destroy(GTK_WIDGET(child->data));
                }
                g_list_free(children);

                // Reload history from the (now empty) file
                load_history(historybar);

                // Force the historybar to update immediately
                gtk_widget_queue_draw(historybar);
            } else {
                std::cerr << "Error: History bar not found." << std::endl;
            }
        }
    } else {
        std::cerr << "Error: widget is not a valid top-level window." << std::endl;
    }
}



///HISTORYBAR
void on_openhistory_button_clicked(GtkWidget *button, GtkWidget *historybar) {
    if (gtk_widget_get_visible(historybar)) {
        gtk_widget_hide(historybar);
    } else {
        gtk_widget_show(historybar);
    }
}
///////
//////END HISTORY

//COOKIES
static void on_delall_cookies_button_clicked(GtkWidget *widget, gpointer data) {
    // Get the top-level window of the widget
    GtkWidget *toplevel = gtk_widget_get_toplevel(widget);

    if (GTK_IS_WIDGET(toplevel) && GTK_IS_WINDOW(toplevel)) {
        // Create a confirmation dialog
        GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(toplevel),
                                                   GTK_DIALOG_MODAL,
                                                   GTK_MESSAGE_QUESTION,
                                                   GTK_BUTTONS_NONE,
                                                   "Are you sure?");
        gtk_widget_set_name(dialog, "dialog");

        // Custom header bar for the dialog
        GtkWidget *header_bar = gtk_header_bar_new();
        gtk_header_bar_set_show_close_button(GTK_HEADER_BAR(header_bar), TRUE);
        gtk_header_bar_set_title(GTK_HEADER_BAR(header_bar), "Delete cookies");
        gtk_widget_set_name(header_bar, "dialog_header_bar");
        gtk_window_set_titlebar(GTK_WINDOW(dialog), header_bar);

        // Add dialog buttons
        GtkWidget *yes_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Yes", GTK_RESPONSE_ACCEPT);
        gtk_widget_set_name(yes_button, "yes_button");
        GtkWidget *no_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "No", GTK_RESPONSE_REJECT);
        gtk_widget_set_name(no_button, "no_button");

        gtk_widget_show_all(dialog);
        gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_CANCEL);

        // Handle the dialog response
        gint response = gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);

        if (response == GTK_RESPONSE_ACCEPT) {
            // Clear the cookies file
            FILE *file = fopen("/usr/share/debris/main_file/cookies.txt", "w");
            if (file != NULL) {
                fclose(file);
            }
        }
        // If "No" is pressed, the dialog is already destroyed, so no further action is needed.
    }
}
//END COOKIES

//CACHES
static void on_delall_caches_button_clicked(GtkWidget *widget, gpointer data) {
    // Get the top-level window of the widget
    GtkWidget *toplevel = gtk_widget_get_toplevel(widget);

    if (GTK_IS_WIDGET(toplevel) && GTK_IS_WINDOW(toplevel)) {
        // Create a confirmation dialog
        GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(toplevel),
                                                   GTK_DIALOG_MODAL,
                                                   GTK_MESSAGE_QUESTION,
                                                   GTK_BUTTONS_NONE,
                                                   "Are you sure?");
        gtk_widget_set_name(dialog, "dialog");

        // Custom header bar for the dialog
        GtkWidget *header_bar = gtk_header_bar_new();
        gtk_header_bar_set_show_close_button(GTK_HEADER_BAR(header_bar), TRUE);
        gtk_header_bar_set_title(GTK_HEADER_BAR(header_bar), "Delete caches");
        gtk_widget_set_name(header_bar, "dialog_header_bar");
        gtk_window_set_titlebar(GTK_WINDOW(dialog), header_bar);

        // Add dialog buttons
        GtkWidget *yes_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Yes", GTK_RESPONSE_ACCEPT);
        gtk_widget_set_name(yes_button, "yes_button");
        GtkWidget *no_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "No", GTK_RESPONSE_REJECT);
        gtk_widget_set_name(no_button, "no_button");

        gtk_widget_show_all(dialog);
        gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_CANCEL);

        // Handle the dialog response
        gint response = gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);

        if (response == GTK_RESPONSE_ACCEPT) {
            // Clear the cacehs file
            FILE *file = fopen("/usr/share/debris/main_file/caches.txt", "w");
            if (file != NULL) {
                fclose(file);
            }
        }
        // If "No" is pressed, the dialog is already destroyed, so no further action is needed.
    }
}
//END CACHES

//SEARCH ENGINE
typedef enum {
    SEARCH_ENGINE_DUCKDUCKGO,
    SEARCH_ENGINE_GOOGLE,
    SEARCH_ENGINE_BING,
    SEARCH_ENGINE_WIKI,
    SEARCH_ENGINE_YAHOO,
    SEARCH_ENGINE_WOLFRAMALPHA,
} SearchEngine;

SearchEngine selected_search_engine = SEARCH_ENGINE_DUCKDUCKGO; // Default search engine


typedef struct {
    WebKitWebView *web_view;
    gchar *current_uri;
    gboolean is_file_uri;
} WebViewData;

void show_confirmation_dialog(GtkButton *button, const gchar *message, void (*on_confirm)(GtkButton *, gpointer), gpointer user_data) {
    GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(gtk_widget_get_toplevel(GTK_WIDGET(button))),
                                               GTK_DIALOG_MODAL,
                                               GTK_MESSAGE_QUESTION,
                                               GTK_BUTTONS_NONE,
                                               "%s", message);
    gtk_widget_set_name(dialog, "dialog");
    
GtkWidget *header_bar = gtk_header_bar_new();
gtk_header_bar_set_show_close_button(GTK_HEADER_BAR(header_bar), TRUE);
gtk_header_bar_set_title(GTK_HEADER_BAR(header_bar), "");
gtk_widget_set_name(header_bar, "dialog_header_bar");

// Add the header bar to the dialog's custom title
gtk_window_set_titlebar(GTK_WINDOW(dialog), header_bar);

// Add buttons
GtkWidget *open_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Yes", GTK_RESPONSE_ACCEPT);
gtk_widget_set_name(open_button, "no_button");
GtkWidget *delete_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "No", GTK_RESPONSE_REJECT);
gtk_widget_set_name(delete_button, "yes_button");

gtk_widget_show_all(dialog);

    gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_CANCEL);

    gint response = gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);

    if (response == GTK_RESPONSE_ACCEPT) {
        on_confirm(button, user_data);
    }
}

//SAVE PREFERENCES ENGINE
void save_engine_choice(SearchEngine engine) {
    const char *engine_str = NULL;
    switch (engine) {
        case SEARCH_ENGINE_DUCKDUCKGO: engine_str = "DuckDuckGo"; break;
        case SEARCH_ENGINE_GOOGLE: engine_str = "Google"; break;
        case SEARCH_ENGINE_BING: engine_str = "Bing"; break;
        case SEARCH_ENGINE_WIKI: engine_str = "Wiki"; break;
        case SEARCH_ENGINE_YAHOO: engine_str = "Yahoo"; break;
        case SEARCH_ENGINE_WOLFRAMALPHA: engine_str = "WolframAlpha"; break;
        default: engine_str = "DuckDuckGo";
    }
    save_preference("search_engine", engine_str);
}

void load_saved_engine_choice(SearchEngine *engine) {
    GKeyFile *key_file = g_key_file_new();
    char *config_path = g_build_filename(g_get_user_config_dir(), "debris_config.ini", NULL);

    *engine = SEARCH_ENGINE_DUCKDUCKGO; // Default

    if (g_key_file_load_from_file(key_file, config_path, G_KEY_FILE_NONE, NULL)) {
        char *engine_str = g_key_file_get_string(key_file, "Preferences", "search_engine", NULL);
        if (engine_str) {
            if (g_strcmp0(engine_str, "Google") == 0) *engine = SEARCH_ENGINE_GOOGLE;
            else if (g_strcmp0(engine_str, "Bing") == 0) *engine = SEARCH_ENGINE_BING;
            else if (g_strcmp0(engine_str, "Wiki") == 0) *engine = SEARCH_ENGINE_WIKI;
            else if (g_strcmp0(engine_str, "Yahoo") == 0) *engine = SEARCH_ENGINE_YAHOO;
            else if (g_strcmp0(engine_str, "WolframAlpha") == 0) *engine = SEARCH_ENGINE_WOLFRAMALPHA;

            printf("Loaded search engine: %s\n", engine_str);
            g_free(engine_str);
        }
    }

    g_free(config_path);
    g_key_file_free(key_file);
}


void confirm_duckduckgo_selection(GtkButton *button, gpointer user_data) {
    selected_search_engine = SEARCH_ENGINE_DUCKDUCKGO;
    save_engine_choice(selected_search_engine);
   // g_print("Selected search engine: DuckDuckGo\n");
}

void confirm_google_selection(GtkButton *button, gpointer user_data) {
    selected_search_engine = SEARCH_ENGINE_GOOGLE;
    save_engine_choice(selected_search_engine);
   // g_print("Selected search engine: Google\n");
}

void confirm_bing_selection(GtkButton *button, gpointer user_data) {
    selected_search_engine = SEARCH_ENGINE_BING;
    save_engine_choice(selected_search_engine);
   // g_print("Selected search engine: Bing\n");
}

void confirm_wiki_selection(GtkButton *button, gpointer user_data) {
    selected_search_engine = SEARCH_ENGINE_WIKI;
    save_engine_choice(selected_search_engine);
   // g_print("Selected search engine: Wiki\n");
}

void confirm_yahoo_selection(GtkButton *button, gpointer user_data) {
    selected_search_engine = SEARCH_ENGINE_YAHOO;
    save_engine_choice(selected_search_engine);
   // g_print("Selected search engine: Wiki\n");
}

void confirm_wolframalpha_selection(GtkButton *button, gpointer user_data) {
    selected_search_engine = SEARCH_ENGINE_WOLFRAMALPHA; 
    save_engine_choice(selected_search_engine);
   // g_print("Selected search engine: Wiki\n");
}


void on_settings_search_duck_clicked(GtkButton *button, gpointer user_data) {
    show_confirmation_dialog(button, "Do you want to use DuckDuckGo?", confirm_duckduckgo_selection, user_data);
}

void on_settings_search_google_clicked(GtkButton *button, gpointer user_data) {
    show_confirmation_dialog(button, "Do you want to use Google?", confirm_google_selection, user_data);
}

void on_settings_search_bing_clicked(GtkButton *button, gpointer user_data) {
    show_confirmation_dialog(button, "Do you want to use Bing?", confirm_bing_selection, user_data);
}

void on_settings_search_wiki_clicked(GtkButton *button, gpointer user_data) {
    show_confirmation_dialog(button, "Do you want to use Wikipedia?", confirm_wiki_selection, user_data);
}

void on_settings_search_yahoo_clicked(GtkButton *button, gpointer user_data) {
    show_confirmation_dialog(button, "Do you want to use yahoo?", confirm_yahoo_selection, user_data);
}

void on_settings_search_wolframalpha_clicked(GtkButton *button, gpointer user_data) {
    show_confirmation_dialog(button, "Do you want to use Wolframalpha?", confirm_wolframalpha_selection, user_data);
}

// MainPage entry activate callback function
void on_mainpage_entry_activate(GtkEntry *entry, gpointer user_data) {
    g_return_if_fail(GTK_IS_NOTEBOOK(user_data));

    GtkNotebook *notebook = GTK_NOTEBOOK(user_data);
    const gchar *input_text = gtk_entry_get_text(entry);

    gchar *escaped_text = g_strdup(input_text);
    for (gchar *p = escaped_text; *p; p++) {
        if (*p == ' ') *p = '+';
    }

    gchar *url;
    if (selected_search_engine == SEARCH_ENGINE_DUCKDUCKGO) {
        url = g_strdup_printf("https://duckduckgo.com/?q=%s&t=ftsa&ia=web", escaped_text);
    } else if (selected_search_engine == SEARCH_ENGINE_GOOGLE) { 
        url = g_strdup_printf("https://www.google.com/search?q=%s", escaped_text);
    } else if (selected_search_engine == SEARCH_ENGINE_BING) { 
        url = g_strdup_printf("https://www.bing.com/search?q=%s", escaped_text);
    } else if (selected_search_engine == SEARCH_ENGINE_WIKI) { 
        url = g_strdup_printf("https://en.wikipedia.org/wiki/%s", escaped_text);
    }  else if (selected_search_engine == SEARCH_ENGINE_YAHOO) { 
        url = g_strdup_printf("https://search.yahoo.com/search?p=%s", escaped_text);
    }  else if (selected_search_engine == SEARCH_ENGINE_WOLFRAMALPHA) { 
        url = g_strdup_printf("https://www.wolframalpha.com/input?i=%s", escaped_text);
    }

    g_free(escaped_text);

   // g_print("Loading URL: %s\n", url);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, url);
    enable_experimental_features(web_view);

    GtkWidget *tab_label = gtk_label_new(url);
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);

    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);

    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}



void on_url_entry_activate(GtkEntry *entry, gpointer user_data) {
    // Ensure the user_data is a GtkNotebook
    g_return_if_fail(GTK_IS_NOTEBOOK(user_data));

    GtkNotebook *notebook = GTK_NOTEBOOK(user_data);
    const gchar *url = gtk_entry_get_text(entry);

   // g_print("Loading URL: %s\n", url);

    // Create a new WebView and load the URL
    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, url);
    enable_experimental_features(web_view);

    // Create a label for the tab
    GtkWidget *tab_label = gtk_label_new(url);

    // Add the new web view to the notebook with the tab label
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);

    // Connect the title-changed signal to update the tab label
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);

    // Show all widgets in the new tab
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}


///OPENFILEINWEBVIEW
static void on_openfilebutton_clicked(GtkWidget *widget, gpointer data) {
    if (!GTK_IS_NOTEBOOK(data)) {
        //g_print("Invalid data parameter: not a GtkNotebook.\n");
        return;
    }

    GtkWidget *open_dialog;
    GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;

    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    //g_print("Widget: %p\n", widget);

    GtkWidget *toplevel = gtk_widget_get_toplevel(widget);

    if (GTK_IS_WIDGET(toplevel) && GTK_IS_WINDOW(toplevel)) {
        //g_print("Widget is a valid top-level window\n");

        open_dialog = gtk_file_chooser_dialog_new("Open File",
                                                  GTK_WINDOW(toplevel),
                                                  action,
                                                  "_Cancel",
                                                  GTK_RESPONSE_CANCEL,
                                                  "_Open",
                                                  GTK_RESPONSE_ACCEPT,
                                                  NULL);
        gtk_widget_set_name(open_dialog, "file_chooser");
        gint res = gtk_dialog_run(GTK_DIALOG(open_dialog));

        if (res == GTK_RESPONSE_ACCEPT) {
            char *filename;
            GtkFileChooser *chooser = GTK_FILE_CHOOSER(open_dialog);
            filename = gtk_file_chooser_get_filename(chooser);

            // Create a new web view and set web view settings
            WebKitWebView *web_view = setup_webview_storage();
            WebKitSettings *settings = webkit_web_view_get_settings(web_view);
            enable_experimental_features(web_view);        
            g_object_set(settings, "allow-file-access-from-file-urls", TRUE, NULL);
            g_object_set(settings, "allow-universal-access-from-file-urls", TRUE, NULL);

           

            // Extract the file extension
            gchar *ext = g_strrstr(filename, ".");

            // Load the appropriate content based on the file extension
            if (ext && g_strcmp0(ext, ".html") == 0) {
                FILE *file = fopen(filename, "r");
                if (file) {
                    fseek(file, 0, SEEK_END);
                    long length = ftell(file);
                    fseek(file, 0, SEEK_SET);
                    char *content = (char *)malloc(length + 1);
                    fread(content, 1, length, file);
                    fclose(file);
                    content[length] = '\0';

                    gchar *dir = g_path_get_dirname(filename);
                    gchar *base_uri = g_filename_to_uri(dir, NULL, NULL);

                    webkit_web_view_load_html(web_view, content, base_uri);

                    free(content);
                    g_free(dir);
                    g_free(base_uri);
                } else {
                    //g_print("Failed to open file.\n");
                }
            } else {
                gchar *file_uri = g_filename_to_uri(filename, NULL, NULL);
                webkit_web_view_load_uri(web_view, file_uri);
                g_free(file_uri);
            }

            WebViewData *web_view_data = g_new0(WebViewData, 1);
            web_view_data->web_view = web_view;
            web_view_data->current_uri = g_strdup(filename);
            web_view_data->is_file_uri = TRUE;

            GtkWidget *tab_label = gtk_label_new(filename);
            int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
            g_object_set_data(G_OBJECT(gtk_notebook_get_nth_page(notebook, page_num)), "web-view-data", web_view_data);

            g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
            g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);

            gtk_notebook_set_current_page(notebook, page_num);

            g_free(filename);
            gtk_widget_show_all(GTK_WIDGET(web_view));
            gtk_widget_show_all(GTK_WIDGET(tab_label));
        }

        gtk_widget_destroy(open_dialog);
    } else {
        //g_print("Error: widget is not a valid top-level window\n");
    }
}

//RELOAD
void on_reload_button_clicked(GtkButton *button, gpointer user_data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(user_data);
    if (GTK_IS_NOTEBOOK(notebook)) {
        int page_num = gtk_notebook_get_current_page(notebook);
        GtkWidget *current_page = gtk_notebook_get_nth_page(notebook, page_num);
        if (WEBKIT_IS_WEB_VIEW(current_page)) {
            WebKitWebView *web_view = WEBKIT_WEB_VIEW(current_page);
            enable_experimental_features(web_view);
            webkit_web_view_reload(web_view);
        } else {
           // g_print("Current page is not a WebKitWebView.\n");
        }
    } else {
      //  g_print("User data is not a GtkNotebook.\n");
    }
}




void on_switch_page(GtkNotebook *notebook, GtkWidget *page, guint page_num, gpointer data) {
    if (WEBKIT_IS_WEB_VIEW(page)) {
        WebKitWebView *web_view = WEBKIT_WEB_VIEW(page);
        const gchar *uri = webkit_web_view_get_uri(web_view);
        GtkWidget *entry = GTK_WIDGET(data);
        gtk_entry_set_text(GTK_ENTRY(entry), uri ? uri : "");
    } else {
        GtkWidget *entry = GTK_WIDGET(data);
        gtk_entry_set_text(GTK_ENTRY(entry), "Main Page");
    }
}

/////CLOSE BUTTON
static void on_close_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(user_data),
                                               GTK_DIALOG_MODAL,
                                               GTK_MESSAGE_QUESTION,
                                               GTK_BUTTONS_NONE, // Use GTK_BUTTONS_NONE to avoid default buttons
                                               "Are you sure?");
    gtk_widget_set_name(dialog, "dialog");

// Add buttons
GtkWidget *open_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Yes", GTK_RESPONSE_YES);
gtk_widget_set_name(open_button, "no_button");
GtkWidget *delete_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "No", GTK_RESPONSE_NO);
gtk_widget_set_name(delete_button, "yes_button");

gtk_widget_show_all(dialog);
    
    // Set default response
    gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_NO);
    
    gint response = gtk_dialog_run(GTK_DIALOG(dialog));
    if (response == GTK_RESPONSE_YES) {
        gtk_widget_destroy(GTK_WIDGET(user_data)); // Close the window
    }
    
    gtk_widget_destroy(dialog); // Destroy the dialog
}
///////



static void on_minimize_button_clicked(GtkWidget *widget, gpointer data) {
    gtk_window_iconify(GTK_WINDOW(data));
}

static void on_maximize_button_clicked(GtkWidget *widget, gpointer data) {
    gtk_window_maximize(GTK_WINDOW(data));
}

static void on_search_button_clicked(GtkWidget *search_button, gpointer user_data) {
    GtkWidget *entry = GTK_WIDGET(user_data);
    
    // Toggle the visibility of the address bar
    if (gtk_widget_get_visible(entry)) {
        gtk_widget_hide(entry); // Hide the address bar
    } else {
        gtk_widget_show(entry); // Show the address bar
        gtk_widget_grab_focus(entry); // Focus on the address bar
        //gtk_button_set_label(GTK_BUTTON(button), ""); // Change button label
    }
}

//CLOSE TAB
void close_active_tab(GtkWidget *widget, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);
 gint page_num = gtk_notebook_get_current_page(notebook);
    gint num_pages = gtk_notebook_get_n_pages(notebook);

    if (page_num != -1 && num_pages > 1) {
        gtk_notebook_remove_page(notebook, page_num);
    }
}


//MENUBAR
void on_menu_button_clicked(GtkWidget *button, GtkWidget *menubar) {
    // Check if the menubar is currently visible
    if (gtk_widget_get_visible(menubar)) {
        // Hide the menubar
        gtk_widget_hide(menubar);
    } else {
        // Show the menubar
        gtk_widget_show(menubar);
    }
}
//////////

//SEARCHENGINEBAR
void on_searchengine_button_clicked(GtkWidget *button, GtkWidget *searchenginebar) {
    if (gtk_widget_get_visible(searchenginebar)) {
        gtk_widget_hide(searchenginebar);
    } else {
        gtk_widget_show(searchenginebar);
    }
}
//////////

//INSTRUMENTBAR
void on_tools_button_clicked(GtkWidget *button, GtkWidget *instrumentbar) {
    if (gtk_widget_get_visible(instrumentbar)) {
        gtk_widget_hide(instrumentbar);
    } else {
        gtk_widget_show(instrumentbar);
    }
}
//////////


//NOTEBAR
void on_note_button_clicked(GtkWidget *button, GtkWidget *notebar) {
    if (gtk_widget_get_visible(notebar)) {
        gtk_widget_hide(notebar);
    } else {
        gtk_widget_show(notebar);
    }
}
//////////

//BOOKMARKBAR
// Function to set notebook data
void set_notebook_data(GtkWidget *bookmarkbar, GtkNotebook *notebook) {
    g_object_set_data(G_OBJECT(bookmarkbar), "notebook", notebook);
}


void on_openbookmark_button_clicked(GtkWidget *button, GtkWidget *bookmarkbar) {
    if (gtk_widget_get_visible(bookmarkbar)) {
        gtk_widget_hide(bookmarkbar);
    } else {
        gtk_widget_show(bookmarkbar);
    }
}

// Global bookmark list
std::vector<std::string> bookmarks;

// Get the correct user-specific bookmark path
std::string get_bookmark_file_path() {
    const gchar *config_dir = g_get_user_config_dir(); // Detect user-specific config directory
    std::string debris_dir = std::string(config_dir) + "/debris";

    // Ensure the ~/.config/debris/ directory exists
    g_mkdir_with_parents(debris_dir.c_str(), 0755);

    return debris_dir + "/bookmarks.txt";
}

// Function to trim the trailing slash from the URL
std::string trim_url(const std::string& url) {
    if (!url.empty() && url.back() == '/') {
        return url.substr(0, url.size() - 1);
    }
    return url;
}

// Function to generate display-friendly label
std::string generate_display_label(const std::string &url) {
    std::string label = url;

    // Remove http:// or https://
    if (label.find("http://") == 0) {
        label.erase(0, 7);
    } else if (label.find("https://") == 0) {
        label.erase(0, 8);
    }

    // Remove www.
    if (label.find("www.") == 0) {
        label.erase(0, 4);
    }

    // Trim text after the domain extension
    std::size_t pos = label.find('/');
    if (pos != std::string::npos) {
        label = label.substr(0, pos); // Keep only the domain
    }

    return label;
}

// Store the full URL as data on the button
void set_full_url(GtkWidget *button, const gchar *url) {
    g_object_set_data_full(G_OBJECT(button), "full-url", g_strdup(url), g_free);
}

const gchar *get_full_url(GtkWidget *button) {
    return (const gchar *)g_object_get_data(G_OBJECT(button), "full-url");
}

// Function to save bookmarks to a user-specific file
void save_bookmarks() {
    std::string bookmark_file = get_bookmark_file_path();
    std::ofstream file(bookmark_file);

    if (!file) {
        std::cerr << "Error opening file for writing: " << bookmark_file << std::endl;
        return;
    }

    for (const std::string &bookmark : bookmarks) {
        file << bookmark << std::endl;
    }

    std::cout << "Bookmarks saved successfully to: " << bookmark_file << std::endl;
}

// Function to remove a bookmark from the user's bookmarks file
void remove_bookmark_from_file(const gchar *url) {
    std::string bookmark_file = get_bookmark_file_path();
    std::string temp_file = bookmark_file + ".temp";

    FILE *file = fopen(bookmark_file.c_str(), "r");
    FILE *temp = fopen(temp_file.c_str(), "w");

    if (!file || !temp) {
        std::cerr << "Error opening bookmark files." << std::endl;
        return;
    }

    gchar line[256];
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, url) == NULL) {
            fputs(line, temp);
        }
    }

    fclose(file);
    fclose(temp);
    remove(bookmark_file.c_str());
    rename(temp_file.c_str(), bookmark_file.c_str());

    std::cout << "Bookmark removed successfully." << std::endl;
}


// Modified callback function to handle the creation of a new page with the bookmark URL or deletion
void on_bookmark_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);
    const gchar *url = get_full_url(button);

    //g_print("Widget: %p\n", button);

    GtkWidget *toplevel = gtk_widget_get_toplevel(button);

    if (GTK_IS_WIDGET(toplevel) && GTK_IS_WINDOW(toplevel)) {
        //g_print("Widget is a valid top-level window\n");

        GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(toplevel),
                                                   GTK_DIALOG_MODAL,
                                                   GTK_MESSAGE_QUESTION,
                                                   GTK_BUTTONS_NONE,
                                                   "You want to open or delete?");
        gtk_widget_set_name(dialog, "dialog");

GtkWidget *header_bar = gtk_header_bar_new();
gtk_header_bar_set_show_close_button(GTK_HEADER_BAR(header_bar), TRUE);
gtk_header_bar_set_title(GTK_HEADER_BAR(header_bar), "");
gtk_widget_set_name(header_bar, "dialog_header_bar");

// Add the header bar to the dialog's custom title
gtk_window_set_titlebar(GTK_WINDOW(dialog), header_bar);

// Add buttons
GtkWidget *open_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Open", GTK_RESPONSE_ACCEPT);
gtk_widget_set_name(open_button, "no_button");
GtkWidget *delete_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Delete", GTK_RESPONSE_REJECT);
gtk_widget_set_name(delete_button, "yes_button");

gtk_widget_show_all(dialog);

        gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_CANCEL);

        gint response = gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);

        if (response == GTK_RESPONSE_ACCEPT) {
            // Create a new web view
            WebKitWebView *web_view = setup_webview_storage();
            webkit_web_view_load_uri(web_view, url);
            enable_experimental_features(web_view);

            // Create a label for the tab
            GtkWidget *tab_label = gtk_label_new(url);

            // Add the new web view to the notebook with the tab label
            int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);

            // Connect the title-changed signal to update the tab label
            g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
            g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);

            // Show all widgets in the new tab
            gtk_widget_show_all(GTK_WIDGET(web_view));
            gtk_widget_show_all(GTK_WIDGET(tab_label));
        } else if (response == GTK_RESPONSE_REJECT) {
            // Remove the bookmark entry from the bookmark bar and the file
            gtk_widget_destroy(button);
            remove_bookmark_from_file(url);
        }
    } else {
        //g_print("Error: widget is not a valid top-level window\n");
    }
}


// Function to load bookmarks from a file
void load_bookmarks(GtkWidget *bookmarkbar) {
    std::ifstream file(get_bookmark_file_path());
    if (!file) {
        std::cerr << "No bookmarks found. Creating an empty file." << std::endl;
        std::ofstream new_file(get_bookmark_file_path());
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;  // Skip empty lines
        }

        bookmarks.push_back(line);

        // Generate display-friendly label for the button
        std::string display_label = generate_display_label(line);

        // Create a new button with the display-friendly label
        GtkWidget *bookmark_button = gtk_button_new_with_label(display_label.c_str());
        gtk_widget_set_name(bookmark_button, "bookmark_button");

        // Align label text
        GtkLabel *label = GTK_LABEL(gtk_bin_get_child(GTK_BIN(bookmark_button)));
        if (GTK_IS_LABEL(label)) {
            gtk_widget_set_halign(GTK_WIDGET(label), GTK_ALIGN_START);
        }

        gtk_box_pack_start(GTK_BOX(bookmarkbar), bookmark_button, FALSE, FALSE, 0);
        gtk_widget_show(bookmark_button);

        // Store the full URL as data
        set_full_url(bookmark_button, line.c_str());

        // Retrieve the notebook and connect the button click signal
        GtkNotebook *notebook = GTK_NOTEBOOK(g_object_get_data(G_OBJECT(bookmarkbar), "notebook"));
        g_signal_connect(bookmark_button, "clicked", G_CALLBACK(on_bookmark_clicked), G_OBJECT(notebook));
    }
}


// Callback function to add bookmark
void add_bookmark(GtkWidget *widget, gpointer data) {
    // Get the notebook widget
    GtkNotebook *notebook = GTK_NOTEBOOK(data);
    int current_page = gtk_notebook_get_current_page(notebook);

    // Get the WebKitWebView for the active tab
    WebKitWebView *web_view = WEBKIT_WEB_VIEW(gtk_notebook_get_nth_page(notebook, current_page));

    if (!WEBKIT_IS_WEB_VIEW(web_view)) {
        std::cerr << "Error: Invalid WebKitWebView pointer." << std::endl;
        return;
    }

    const gchar *uri = webkit_web_view_get_uri(web_view);
    if (uri) {
        std::string bookmark = trim_url(uri);

        // Add bookmark to the list
        bookmarks.push_back(bookmark);

        // **Save bookmarks to the user's config folder**
        std::ofstream file(get_bookmark_file_path());
        if (file) {
            for (const std::string &b : bookmarks) {
                file << b << std::endl;
            }
            file.close();
        } else {
            std::cerr << "Error writing bookmarks to file." << std::endl;
        }

        // Generate display-friendly label
        std::string display_label = generate_display_label(bookmark);

        // Create a button with the label
        GtkWidget *bookmark_button = gtk_button_new_with_label(display_label.c_str());
        gtk_widget_set_name(bookmark_button, "bookmark_button");
        
        GtkLabel *label = GTK_LABEL(gtk_bin_get_child(GTK_BIN(bookmark_button)));
        gtk_widget_set_halign(GTK_WIDGET(label), GTK_ALIGN_START); 

        GtkWidget *bookmarkbar = GTK_WIDGET(g_object_get_data(G_OBJECT(widget), "bookmarkbar"));
        if (bookmarkbar) {
            gtk_box_pack_start(GTK_BOX(bookmarkbar), bookmark_button, FALSE, FALSE, 0);
            gtk_widget_show(bookmark_button);

            // Store the full URL as data on the button
            set_full_url(bookmark_button, uri);

            // Retrieve notebook reference and connect the signal
            GtkNotebook *notebook = GTK_NOTEBOOK(g_object_get_data(G_OBJECT(bookmarkbar), "notebook"));
            g_signal_connect(bookmark_button, "clicked", G_CALLBACK(on_bookmark_clicked), notebook);
        } else {
            std::cerr << "Error: Bookmark bar not found." << std::endl;
        }
    } else {
        std::cerr << "Error: No URI found for the current web view." << std::endl;
    }
}


//Delete all bookmarks
static void on_delall_bookmarks_button_clicked(GtkWidget *widget, gpointer data) {
    GtkWidget *toplevel = gtk_widget_get_toplevel(widget);

    if (GTK_IS_WIDGET(toplevel) && GTK_IS_WINDOW(toplevel)) {
        // Confirmation dialog
        GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(toplevel),
                                                   GTK_DIALOG_MODAL,
                                                   GTK_MESSAGE_QUESTION,
                                                   GTK_BUTTONS_NONE,
                                                   "Are you sure?");
        gtk_widget_set_name(dialog, "dialog");

        // Custom header bar
        GtkWidget *header_bar = gtk_header_bar_new();
        gtk_header_bar_set_show_close_button(GTK_HEADER_BAR(header_bar), TRUE);
        gtk_header_bar_set_title(GTK_HEADER_BAR(header_bar), "Delete bookmarks");
        gtk_widget_set_name(header_bar, "dialog_header_bar");
        gtk_window_set_titlebar(GTK_WINDOW(dialog), header_bar);

        // Add buttons
        GtkWidget *yes_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Yes", GTK_RESPONSE_ACCEPT);
        gtk_widget_set_name(yes_button, "yes_button");
        GtkWidget *no_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "No", GTK_RESPONSE_REJECT);
        gtk_widget_set_name(no_button, "no_button");

        gtk_widget_show_all(dialog);
        gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_CANCEL);

        gint response = gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);

        if (response == GTK_RESPONSE_ACCEPT) {
            // Get the user-specific bookmarks file path
            std::string bookmark_file = get_bookmark_file_path();

            // Empty the file instead of deleting
            std::ofstream file(bookmark_file, std::ios::trunc);
            if (file) {
                file.close();
            } else {
                std::cerr << "Error clearing bookmarks file: " << bookmark_file << std::endl;
            }
        }
    }
}



////END BOOKMRKS

//////URL////
static void load_url(WebKitWebView *web_view, const gchar *url) {
    // Load the specified URL in the WebView
    webkit_web_view_load_uri(web_view, url);
}


static void on_google_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.google.com");
    enable_experimental_features(web_view);

    GtkWidget *tab_label = gtk_label_new("Google");
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);

    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}


static void on_bing_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.bing.com");
    enable_experimental_features(web_view);
    
    GtkWidget *tab_label = gtk_label_new("Bing");
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_wiki_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.wikipedia.com");
    enable_experimental_features(web_view);
    
    GtkWidget *tab_label = gtk_label_new("Wikipedia");
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_duck_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://duckduckgo.com");
    GtkWidget *tab_label = gtk_label_new("DuckDuckGo");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_yahoo_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://search.yahoo.com/");
    GtkWidget *tab_label = gtk_label_new("Yahoo");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_wolframalpha_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.wolframalpha.com/");
    GtkWidget *tab_label = gtk_label_new("WolframAlpha");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

//NEWS URL
static void on_ap_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://apnews.com/");
    GtkWidget *tab_label = gtk_label_new("AP NEWS");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_alarabya_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://english.alarabiya.net/");
    GtkWidget *tab_label = gtk_label_new("Al Arabya");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_aljazeera_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.aljazeera.com/");
    GtkWidget *tab_label = gtk_label_new("Aljazeera");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_ecns_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "http://www.ecns.cn/");
    GtkWidget *tab_label = gtk_label_new("ECNS");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_moscowtimes_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.themoscowtimes.com/");
    GtkWidget *tab_label = gtk_label_new("M T");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_nytimes_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.nytimes.com/international/");
    GtkWidget *tab_label = gtk_label_new("NYT");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_cnbc_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.cnbc.com/world/?region=world");
    GtkWidget *tab_label = gtk_label_new("CNBC");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_lemonde_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.lemonde.fr/en/");
    GtkWidget *tab_label = gtk_label_new("Le Monde");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_euronews_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.euronews.com/");
    GtkWidget *tab_label = gtk_label_new("Euronews");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_nature_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.nature.com/");
    GtkWidget *tab_label = gtk_label_new("NATURE");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_phys_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://phys.org/");
    GtkWidget *tab_label = gtk_label_new("PHYS.ORG");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_sciencenews_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.sciencenews.org/");
    GtkWidget *tab_label = gtk_label_new("SCIENCE");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}
// VIDEO URL
static void on_arte_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.arte.tv/en/");
    GtkWidget *tab_label = gtk_label_new("ARTE");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_bbc_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.bbc.com/video");
    GtkWidget *tab_label = gtk_label_new("BBC");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_docplus_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.docplus.com/");
    GtkWidget *tab_label = gtk_label_new("DOC+");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_dailymotion_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.dailymotion.com/");
    GtkWidget *tab_label = gtk_label_new("Dailymotion");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_jove_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.jove.com/");
    GtkWidget *tab_label = gtk_label_new("Jove");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_nasatv_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://plus.nasa.gov/");
    GtkWidget *tab_label = gtk_label_new("Nasa+");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_odysee_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://odysee.com/");
    GtkWidget *tab_label = gtk_label_new("Odysee");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_pbs_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.pbs.org/");
    GtkWidget *tab_label = gtk_label_new("PBS");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_scienceaaas_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.science.org/videos");
    GtkWidget *tab_label = gtk_label_new("Science");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_ted_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.ted.com/");
    GtkWidget *tab_label = gtk_label_new("TED");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_twitch_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.twitch.tv/");
    GtkWidget *tab_label = gtk_label_new("Twitch");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_youtube_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.youtube.com");
    GtkWidget *tab_label = gtk_label_new("YouTube");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}
//start link education
static void on_biodigital_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.biodigital.com/");
    GtkWidget *tab_label = gtk_label_new("Biodigital");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_dlmf_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://dlmf.nist.gov/");
    GtkWidget *tab_label = gtk_label_new("DLMF");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_efunda_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.efunda.com/home.cfm");
    GtkWidget *tab_label = gtk_label_new("EFUNDA");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_mit_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://ocw.mit.edu/");
    GtkWidget *tab_label = gtk_label_new("MIT");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_onezoom_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.onezoom.org/");
    GtkWidget *tab_label = gtk_label_new("OneZoom");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_pbdb_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://paleobiodb.org/#/");
    GtkWidget *tab_label = gtk_label_new("PBDB");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_pt_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://ptable.com/#Properties");
    GtkWidget *tab_label = gtk_label_new("PT");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_semantic_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.semanticscholar.org/");
    GtkWidget *tab_label = gtk_label_new("Semantic");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_sciencedirect_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.sciencedirect.com/");
    GtkWidget *tab_label = gtk_label_new("SDirect");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_sciencegov_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.science.gov/");
    GtkWidget *tab_label = gtk_label_new("S.gov");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_wwh_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://www.worldhistory.org/");
    GtkWidget *tab_label = gtk_label_new("WHE");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}

static void on_wwt_button_clicked(GtkWidget *button, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    WebKitWebView *web_view = setup_webview_storage();
    webkit_web_view_load_uri(web_view, "https://worldwidetelescope.org/webclient/");
    GtkWidget *tab_label = gtk_label_new("WWT");
    enable_experimental_features(web_view);
    
    int page_num = gtk_notebook_append_page(notebook, GTK_WIDGET(web_view), tab_label);
    g_signal_connect(web_view, "notify::title", G_CALLBACK(on_title_changed), GINT_TO_POINTER(page_num));
    g_signal_connect(web_view, "notify::uri", G_CALLBACK(on_uri_changed), notebook);
    gtk_widget_show_all(GTK_WIDGET(web_view));
    gtk_widget_show_all(GTK_WIDGET(tab_label));
}
//END URL

/////NOTEAPP FUNCTION
// NEW TAB
GdkPixbuf *saved_pixbuf = NULL;
static gboolean on_paste_clipboard(GtkWidget *text_view, GdkEvent *event, gpointer data);

void create_new_tab(GtkWidget *widget, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    // Create a scrolled window
    GtkWidget *scrolled_window = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_window), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

    // Create a text view widget
    GtkWidget *text_view = gtk_text_view_new();
    gtk_widget_set_name(text_view, "text_view");

    // Set maximum width to accommodate approximately 48 characters
    int char_width = 12; // Adjust this based on your font
    int max_width = char_width * 48; // 100 characters
    gtk_widget_set_size_request(text_view, max_width, -1); // -1 means no limit on height

    // Set maximum height to 900 pixels
    gtk_widget_set_size_request(scrolled_window, -1, 900); // -1 means no limit on width

    // Enable text wrapping
    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view), GTK_WRAP_WORD);
    // Set some sample text
    GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
    gtk_text_buffer_set_text(buffer, "", -1);
    // Add the text view to the scrolled window
    gtk_container_add(GTK_CONTAINER(scrolled_window), text_view);
    // Create a tab label
    GtkWidget *tab_label = gtk_label_new("New Tab");
    // Add the new page to the notebook
    gtk_notebook_append_page(notebook, scrolled_window, tab_label);
    gtk_widget_show_all(GTK_WIDGET(notebook));

    // Connect the paste-clipboard event
    g_signal_connect(text_view, "paste-clipboard", G_CALLBACK(on_paste_clipboard), NULL);
}

//////////
//CLOSETAB
void close_active_notetab(GtkWidget *widget, gpointer data) {
    GtkNotebook *notebook = GTK_NOTEBOOK(data);
    gint page_num = gtk_notebook_get_current_page(notebook);
    gint num_pages = gtk_notebook_get_n_pages(notebook);

    if (page_num != -1 && num_pages > 1) {
        gtk_notebook_remove_page(notebook, page_num);
    }
}

///////////////

//CENTER TEXT
static void on_center_text(GtkWidget *widget, gpointer data) {
    // Ensure the data parameter is a GtkNotebook
    if (!GTK_IS_NOTEBOOK(data)) {
        //g_print("Invalid data parameter: not a GtkNotebook.\n");
        return;
    }

    // Get the GtkNotebook from the data parameter
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    // Get the current page number
    gint page_num = gtk_notebook_get_current_page(notebook);
    if (page_num == -1) {
        //g_print("No active tab found.\n");
        return;
    }

    // Get the scrolled window in the current page
    GtkWidget *scrolled_window = gtk_notebook_get_nth_page(notebook, page_num);
    if (!GTK_IS_SCROLLED_WINDOW(scrolled_window)) {
        //g_print("Scrolled window not found in the current page.\n");
        return;
    }

    // Get the text view from the scrolled window
    GtkWidget *text_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
    if (!GTK_IS_TEXT_VIEW(text_view)) {
        //g_print("Text view not found in the scrolled window.\n");
        return;
    }

    GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
    GtkTextIter start, end;

    // Get the selected text start and end
    if (gtk_text_buffer_get_selection_bounds(buffer, &start, &end)) {
        // Get the tag table of the text buffer
        GtkTextTagTable *tag_table = gtk_text_buffer_get_tag_table(buffer);

        // Check if the tag already exists
        GtkTextTag *centered_tag = gtk_text_tag_table_lookup(tag_table, "centered");
        if (!centered_tag) {
            // Create the tag only if it doesn't already exist
            centered_tag = gtk_text_buffer_create_tag(buffer, "centered", "justification", GTK_JUSTIFY_CENTER, NULL);
        }

        // Apply the tag to the selected text
        gtk_text_buffer_apply_tag(buffer, centered_tag, &start, &end);
    } else {
        //g_print("No text selected.\n");
    }
}

///////
////ALIGN LEFT
static void on_left_align_text(GtkWidget *widget, gpointer data) {
    // Ensure the data parameter is a GtkNotebook
    if (!GTK_IS_NOTEBOOK(data)) {
        //g_print("Invalid data parameter: not a GtkNotebook.\n");
        return;
    }

    // Get the GtkNotebook from the data parameter
    GtkNotebook *notebook = GTK_NOTEBOOK(data);

    // Get the current page number
    gint page_num = gtk_notebook_get_current_page(notebook);
    if (page_num == -1) {
        //g_print("No active tab found.\n");
        return;
    }

    // Get the scrolled window in the current page
    GtkWidget *scrolled_window = gtk_notebook_get_nth_page(notebook, page_num);
    if (!GTK_IS_SCROLLED_WINDOW(scrolled_window)) {
        //g_print("Scrolled window not found in the current page.\n");
        return;
    }

    // Get the text view from the scrolled window
    GtkWidget *text_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
    if (!GTK_IS_TEXT_VIEW(text_view)) {
        //g_print("Text view not found in the scrolled window.\n");
        return;
    }

    GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
    GtkTextIter start, end;

    // Get the selected text start and end
    if (gtk_text_buffer_get_selection_bounds(buffer, &start, &end)) {
        // Get the tag table of the text buffer
        GtkTextTagTable *tag_table = gtk_text_buffer_get_tag_table(buffer);

        // Check if the tag already exists
        GtkTextTag *left_align_tag = gtk_text_tag_table_lookup(tag_table, "left_aligned");
        if (!left_align_tag) {
            // Create the tag only if it doesn't already exist
            left_align_tag = gtk_text_buffer_create_tag(buffer, "left_aligned", "justification", GTK_JUSTIFY_LEFT, NULL);
        }

        // Apply the tag to the selected text
        gtk_text_buffer_apply_tag(buffer, left_align_tag, &start, &end);
    } else {
        //g_print("No text selected.\n");
    }
}

//////



//LIST
static gboolean in_list_mode = FALSE;
static gboolean last_enter = FALSE;

static void on_list_text(GtkWidget *widget, gpointer data) {
    if (!GTK_IS_NOTEBOOK(data)) {
        return;
    }

    GtkNotebook *notebook = GTK_NOTEBOOK(data);
    gint page_num = gtk_notebook_get_current_page(notebook);
    if (page_num == -1) {
        return;
    }

    GtkWidget *scrolled_window = gtk_notebook_get_nth_page(notebook, page_num);
    if (!GTK_IS_SCROLLED_WINDOW(scrolled_window)) {
        return;
    }

    GtkWidget *text_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
    if (!GTK_IS_TEXT_VIEW(text_view)) {
        return;
    }

    GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
    GtkTextIter iter;
    gtk_text_buffer_get_end_iter(buffer, &iter);

    // Insert bullet point
    gtk_text_buffer_insert(buffer, &iter, "• ", -1);

    // Move cursor to the correct position after the bullet
    gtk_text_buffer_get_end_iter(buffer, &iter);
    gtk_text_buffer_place_cursor(buffer, &iter);

    // Ensure focus is set on the text view to allow immediate typing
    gtk_widget_grab_focus(GTK_WIDGET(text_view));

    in_list_mode = TRUE;
    last_enter = FALSE;
}



static gboolean on_key_press(GtkWidget *widget, GdkEventKey *event, gpointer data) {
    if (event->keyval != GDK_KEY_Return) {
        last_enter = FALSE;
        return FALSE; // Allow normal key behavior
    }

    GtkNotebook *notebook = GTK_NOTEBOOK(data);
    gint page_num = gtk_notebook_get_current_page(notebook);
    if (page_num == -1) {
        return FALSE;
    }

    GtkWidget *scrolled_window = gtk_notebook_get_nth_page(notebook, page_num);
    if (!GTK_IS_SCROLLED_WINDOW(scrolled_window)) {
        return FALSE;
    }

    GtkWidget *text_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
    if (!GTK_IS_TEXT_VIEW(text_view)) {
        return FALSE;
    }

    GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
    GtkTextIter iter;
    gtk_text_buffer_get_end_iter(buffer, &iter);

    // **Ensure Enter adds a new bullet every time**
    if (in_list_mode) {
        gtk_text_buffer_insert(buffer, &iter, "\n• ", -1);
        gtk_text_buffer_get_end_iter(buffer, &iter);
        gtk_text_view_place_cursor_onscreen(GTK_TEXT_VIEW(text_view));
        return TRUE; // Prevent default Enter behavior
    } else {
        // If not in list mode, just add a newline
        gtk_text_buffer_insert(buffer, &iter, "\n", -1);
    }

    return TRUE; // Consume Enter key event
}

static gboolean on_enter_pressed(GtkWidget *widget, GdkEventKey *event, gpointer data) {
    if (event->keyval != GDK_KEY_Return) {
        last_enter = FALSE;
        return FALSE; // Allow normal key behavior
    }

    if (!GTK_IS_NOTEBOOK(data)) {
        return FALSE;
    }

    GtkNotebook *notebook = GTK_NOTEBOOK(data);
    gint page_num = gtk_notebook_get_current_page(notebook);
    if (page_num == -1) {
        return FALSE;
    }

    GtkWidget *scrolled_window = gtk_notebook_get_nth_page(notebook, page_num);
    if (!GTK_IS_SCROLLED_WINDOW(scrolled_window)) {
        return FALSE;
    }

    GtkWidget *text_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
    if (!GTK_IS_TEXT_VIEW(text_view)) {
        return FALSE;
    }

    GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
    GtkTextIter iter;

    // Get the current cursor position
    GtkTextMark *insert_mark = gtk_text_buffer_get_insert(buffer);
    gtk_text_buffer_get_iter_at_mark(buffer, &iter, insert_mark);

    // Check if the current line contains a bullet (for list mode behavior)
    GtkTextIter start_iter;
    gtk_text_buffer_get_iter_at_line(buffer, &start_iter, gtk_text_iter_get_line(&iter));
    gchar *line_text = gtk_text_buffer_get_text(buffer, &start_iter, &iter, FALSE);

    if (in_list_mode) {
        if (last_enter) {
            // If Enter is pressed twice without typing, remove the last bullet instead of inserting a new line
            gtk_text_buffer_delete(buffer, &start_iter, &iter);
            in_list_mode = FALSE;
        } else {
            // Insert a new bullet point at the cursor position
            gtk_text_buffer_insert(buffer, &iter, "\n• ", -1);
            last_enter = TRUE;
        }
    } else {
        // Insert a new line at the cursor position and move text below
        gtk_text_buffer_insert(buffer, &iter, "\n", -1);
    }

    g_free(line_text);
    gtk_text_view_place_cursor_onscreen(GTK_TEXT_VIEW(text_view));
    gtk_widget_grab_focus(GTK_WIDGET(text_view));

    return TRUE; // Prevent default Enter behavior
}




////// END LIST

///////OPENFILE////////


char** read_words_from_file(const char *filename, int *num_words) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Unable to open file");
        return NULL;
    }

    char **words = NULL;
    char buffer[100];
    *num_words = 0;

    while (fgets(buffer, sizeof(buffer), file)) {
        buffer[strcspn(buffer, "\n")] = 0;  // Remove newline character
        words = (char **)realloc(words, sizeof(char *) * (*num_words + 1));
        words[*num_words] = strdup(buffer);
        (*num_words)++;
    }

    fclose(file);
    return words;
}
// Function to read keywords from a file and store them in a GList
GList *read_keywords_from_file(const char *filename) {
    int num_words;
    char **words = read_words_from_file(filename, &num_words);
    if (!words) return NULL;

    GList *keywords = NULL;
    for (int i = 0; i < num_words; i++) {
        keywords = g_list_append(keywords, words[i]);
    }
    return keywords;
}

// Function to highlight words in the buffer
void highlight_words_in_buffer(GtkTextBuffer *buffer, GtkTextTag *tag, const char **words, int num_words) {
    GtkTextIter start, match_start, match_end;

    for (int i = 0; i < num_words; i++) {
        const char *word = words[i];
        gtk_text_buffer_get_start_iter(buffer, &start);
        while (gtk_text_iter_forward_search(&start, word, GTK_TEXT_SEARCH_TEXT_ONLY, &match_start, &match_end, NULL)) {
            gtk_text_buffer_apply_tag(buffer, tag, &match_start, &match_end);
            start = match_end; // Move start to the end of the last match
        }
    }
}
// Function to highlight keywords in the buffer
void highlight_keywords(GtkTextBuffer *buffer, GList *keywords) {
    // Create the highlight tags
    GtkTextTag *tag_orange = gtk_text_buffer_create_tag(buffer, "highlight_orange", "foreground", "#ff7800", NULL);
    GtkTextTag *tag_red = gtk_text_buffer_create_tag(buffer, "highlight_red", "foreground", "#33b2a4", NULL);
    GtkTextTag *tag_orange2 = gtk_text_buffer_create_tag(buffer, "highlight_orange2", "foreground", "#FB8354", NULL);
    GtkTextTag *tag_blue = gtk_text_buffer_create_tag(buffer, "highlight_blue", "foreground", "#3CA1F8", NULL);    
    GtkTextTag *tag_green = gtk_text_buffer_create_tag(buffer, "highlight_green", "foreground", "#33b2a4", NULL);   
 
    gboolean is_css = FALSE;
    gboolean is_symbol = FALSE;
    gboolean is_blue = FALSE;
    gboolean is_library = FALSE;

    for (GList *l = keywords; l != NULL; l = l->next) {
        const char *word = (const char *)l->data;

        if (g_strcmp0(word, "CSS") == 0) {
            is_css = TRUE;
            is_symbol = FALSE;
            is_blue = FALSE;
            is_library = FALSE;
        } else if (g_strcmp0(word, "SYMBOL") == 0) {
            is_symbol = TRUE;
            is_css = FALSE;
            is_blue = FALSE;
            is_library = FALSE;
        } else if (g_strcmp0(word, "BLUE") == 0) {
            is_symbol = FALSE;
            is_css = FALSE;
            is_blue = TRUE;
            is_library = FALSE;
        } else if (g_strcmp0(word, "LIBRARY") == 0) {
            is_symbol = FALSE;
            is_css = FALSE;
            is_blue = FALSE;
            is_library = TRUE;
        }
        
        

        if (is_css) {
            highlight_words_in_buffer(buffer, tag_red, &word, 1);
        } else if (is_symbol) {
            highlight_words_in_buffer(buffer, tag_orange2, &word, 1);
        } else if (is_blue) {
            highlight_words_in_buffer(buffer, tag_blue, &word, 1);
        } else if (is_library) {
            highlight_words_in_buffer(buffer, tag_green, &word, 1);
        } else {
            highlight_words_in_buffer(buffer, tag_orange, &word, 1);
        }
    }
}

// Declare Boolean format flags
gboolean format_css = FALSE;
gboolean format_html = FALSE;
gboolean format_c = FALSE;
gboolean format_js = FALSE;

static void on_buttonopentext_clicked(GtkWidget *widget, gpointer data) {
    if (!GTK_IS_NOTEBOOK(data)) {
        return;
    }

    GtkWidget *open_dialog_note;
    GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
    GtkNotebook *notebook = GTK_NOTEBOOK(data);
    gint page_num = gtk_notebook_get_current_page(notebook);
    if (page_num == -1) {
        return;
    }

    GtkWidget *scrolled_window = gtk_notebook_get_nth_page(notebook, page_num);
    if (!GTK_IS_SCROLLED_WINDOW(scrolled_window)) {
        return;
    }

    GtkWidget *text_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
    if (!GTK_IS_TEXT_VIEW(text_view)) {
        return;
    }

    open_dialog_note = gtk_file_chooser_dialog_new("Open File",
                                                   GTK_WINDOW(gtk_widget_get_toplevel(widget)),
                                                   action, "_Cancel", GTK_RESPONSE_CANCEL,
                                                   "_Open", GTK_RESPONSE_ACCEPT, NULL);
    gtk_widget_set_name(open_dialog_note, "file_chooser");
    gint res = gtk_dialog_run(GTK_DIALOG(open_dialog_note));

    if (res == GTK_RESPONSE_ACCEPT) {
        char *filename;
        GtkFileChooser *chooser = GTK_FILE_CHOOSER(open_dialog_note);
        filename = gtk_file_chooser_get_filename(chooser);

        char tab_label[11]; 
        strncpy(tab_label, g_path_get_basename(filename), 10);
        tab_label[10] = '\0';

        GtkWidget *label = gtk_label_new(tab_label);
        gtk_notebook_set_tab_label(notebook, scrolled_window, label);

        FILE *file = fopen(filename, "r");
        if (file) {
            fseek(file, 0, SEEK_END);
            long length = ftell(file);
            fseek(file, 0, SEEK_SET);
            char *content = (char *)malloc(length + 1);
            fread(content, 1, length, file);
            fclose(file);
            content[length] = '\0';

            GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
            gtk_text_buffer_set_text(buffer, content, -1);

            const char *ext = strrchr(filename, '.');
            if (ext) {
                GList *keywords = NULL;

                format_css = format_html = format_c = format_js = FALSE; // Reset format flags

                if (g_strcmp0(ext, ".html") == 0 || g_strcmp0(ext, ".htm") == 0) {
                    keywords = read_keywords_from_file("/usr/share/debris/main_file/keywordsht.txt");
                    format_html = TRUE;
                } else if (g_strcmp0(ext, ".c") == 0 || g_strcmp0(ext, ".cpp") == 0) {
                    keywords = read_keywords_from_file("/usr/share/debris/main_file/keywordsc.txt");
                    format_c = TRUE;
                } else if (g_strcmp0(ext, ".js") == 0 || g_strcmp0(ext, ".json") == 0) {
                    keywords = read_keywords_from_file("/usr/share/debris/main_file/keywordsj.txt");
                    format_js = TRUE;
                } else if (g_strcmp0(ext, ".css") == 0) {
                    keywords = read_keywords_from_file("/usr/share/debris/main_file/keywordscss.txt");
                    format_css = TRUE;
                }

                if (keywords) {
                    highlight_keywords(buffer, keywords);
                    g_list_free_full(keywords, g_free);
                }

                g_print("Detected format: %s\n",
                        format_html ? "html" :
                        format_css ? "css" :
                        format_c ? "c" :
                        format_js ? "js" : "Unknown");
            }

            g_free(content);
            g_free(filename);
            gtk_widget_destroy(open_dialog_note);
        }
    }
}

// SAVE AS FILE
// Function to save notes (text + images) as a PDF

static gboolean image_predicate(gunichar ch, gpointer user_data) {
    return ch == 0xFFFC; // Detect image placeholders
}

static void on_buttonsavetext_clicked(GtkWidget *widget, gpointer data) {
    if (!GTK_IS_NOTEBOOK(data)) {
        return;
    }

    GtkWidget *toplevel = gtk_widget_get_toplevel(widget);
    if (!GTK_IS_WIDGET(toplevel) || !GTK_IS_WINDOW(toplevel)) {
        return;
    }

    GtkWidget *dialog;
    GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_SAVE;
    GtkNotebook *notebook = GTK_NOTEBOOK(data);
    gint page_num = gtk_notebook_get_current_page(notebook);
    if (page_num == -1) {
        return;
    }

    GtkWidget *scrolled_window = gtk_notebook_get_nth_page(notebook, page_num);
    if (!GTK_IS_SCROLLED_WINDOW(scrolled_window)) {
        return;
    }

    GtkWidget *text_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
    if (!GTK_IS_TEXT_VIEW(text_view)) {
        return;
    }

    dialog = gtk_file_chooser_dialog_new("Save ODT File",
                                         GTK_WINDOW(toplevel),
                                         action,
                                         "_Cancel",
                                         GTK_RESPONSE_CANCEL,
                                         "_Save",
                                         GTK_RESPONSE_ACCEPT,
                                         NULL);
    gtk_widget_set_name(dialog, "file_chooser");
    gint res = gtk_dialog_run(GTK_DIALOG(dialog));

    if (res == GTK_RESPONSE_ACCEPT) {
        char *filename;
        GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
        filename = gtk_file_chooser_get_filename(chooser);

        GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
        GtkTextIter start, end;
        gtk_text_buffer_get_start_iter(buffer, &start);
        gtk_text_buffer_get_end_iter(buffer, &end);

        // Open file for writing ODT content
        FILE *file = fopen(filename, "w");
        if (!file) {
            g_free(filename);
            gtk_widget_destroy(dialog);
            return;
        }

        // Iterate through text and handle images inline
        GtkTextIter iter = start;
        while (!gtk_text_iter_is_end(&iter)) {
            gunichar ch = gtk_text_iter_get_char(&iter);

            if (ch == 0xFFFC) { // Detect images
                GdkPixbuf *pixbuf = gtk_text_iter_get_pixbuf(&iter);
                if (pixbuf) {
                    // Set image path to match the ODT file directory
char image_filename[256];
sprintf(image_filename, "%s/image_%ld.png", g_path_get_dirname(filename), g_get_real_time());
gdk_pixbuf_save(pixbuf, image_filename, "png", NULL, NULL);


                    // Replace placeholder with actual image reference
                    fprintf(file, "[Image: %s]", image_filename);
                }
            } else {
                gchar char_str[6]; // UTF-8 characters can be up to 6 bytes
                g_unichar_to_utf8(ch, char_str);
                fprintf(file, "%s", char_str);
            }

            gtk_text_iter_forward_char(&iter);
        }

        fprintf(file, "\n"); // Ensure proper line break
        fclose(file);

        // Update notebook tab label
        char tab_label[11]; // 10 characters + null terminator
        strncpy(tab_label, g_path_get_basename(filename), 10);
        tab_label[10] = '\0';
        GtkWidget *label = gtk_label_new(tab_label);
        gtk_notebook_set_tab_label(notebook, scrolled_window, label);

        // Apply syntax highlighting
        const char *ext = strrchr(filename, '.');
        if (ext) {
            GList *keywords = NULL;

            format_html = format_css = format_c = format_js = FALSE; // Reset flags

            if (g_strcmp0(ext, ".html") == 0 || g_strcmp0(ext, ".htm") == 0) {
                keywords = read_keywords_from_file("/usr/share/debris/main_file/keywordsht.txt");
                format_html = TRUE;
            } else if (g_strcmp0(ext, ".c") == 0 || g_strcmp0(ext, ".cpp") == 0) {
                keywords = read_keywords_from_file("/usr/share/debris/main_file/keywordsc.txt");
                format_c = TRUE;
            } else if (g_strcmp0(ext, ".js") == 0 || g_strcmp0(ext, ".json") == 0) {
                keywords = read_keywords_from_file("/usr/share/debris/main_file/keywordsj.txt");
                format_js = TRUE;
            } else if (g_strcmp0(ext, ".css") == 0) {
                keywords = read_keywords_from_file("/usr/share/debris/main_file/keywordscss.txt");
                format_css = TRUE;
            }

            if (keywords) {
                highlight_keywords(buffer, keywords);
                g_list_free_full(keywords, g_free);
            }
        }

        g_print("Detected format: %s\n",
                format_html ? "html" :
                format_css ? "css" :
                format_c ? "c" :
                format_js ? "js" : "Unknown");

        g_free(filename);
    }

    gtk_widget_destroy(dialog);
}


//END SAVE AS

//COMMENT
static void on_comment_text(GtkWidget *widget, gpointer data) {
    if (!GTK_IS_NOTEBOOK(data)) {
        return;
    }

    GtkNotebook *notebook = GTK_NOTEBOOK(data);
    gint page_num = gtk_notebook_get_current_page(notebook);
    if (page_num == -1) {
        return;
    }

    GtkWidget *scrolled_window = gtk_notebook_get_nth_page(notebook, page_num);
    if (!GTK_IS_SCROLLED_WINDOW(scrolled_window)) {
        return;
    }

    GtkWidget *text_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
    if (!GTK_IS_TEXT_VIEW(text_view)) {
        return;
    }

    GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
    GtkTextIter iter;

    // Get cursor position correctly
    GtkTextMark *insert_mark = gtk_text_buffer_get_insert(buffer);
    gtk_text_buffer_get_iter_at_mark(buffer, &iter, insert_mark);

    // Ensure a valid format is detected
    if (!format_html && !format_css && !format_c && !format_js) {
        g_print("No format detected\n");
        return;
    }

    g_print("Detected format: %s\n",
            format_html ? "html" :
            format_css ? "css" :
            format_c ? "c" :
            format_js ? "js" : "Unknown");

    const char *comment_syntax = "/* Comment */"; // Default to C-style comments

    // If format is HTML, check if inside <style> or <script>
    if (format_html) {
        GtkTextIter start_iter, end_iter;
        gtk_text_buffer_get_start_iter(buffer, &start_iter);
        gtk_text_buffer_get_end_iter(buffer, &end_iter);

        gchar *full_text = gtk_text_buffer_get_text(buffer, &start_iter, &end_iter, FALSE);
        gint cursor_offset = gtk_text_iter_get_offset(&iter);

        gboolean in_style_tag = FALSE, in_script_tag = FALSE;
        if (full_text) {
            gchar *style_pos = g_strrstr_len(full_text, cursor_offset, "<style>");
            gchar *script_pos = g_strrstr_len(full_text, cursor_offset, "<script>");
            gchar *style_close = g_strrstr_len(full_text, cursor_offset, "</style>");
            gchar *script_close = g_strrstr_len(full_text, cursor_offset, "</script>");

            if (style_pos && (!style_close || style_pos > style_close)) {
                in_style_tag = TRUE;
            }
            if (script_pos && (!script_close || script_pos > script_close)) {
                in_script_tag = TRUE;
            }

            g_free(full_text);
        }

        // Use <!-- Comment --> if NOT inside <style> or <script>
        if (!in_style_tag && !in_script_tag) {
            comment_syntax = "<!-- Comment -->";
        }
    }

    g_print("Comment syntax: %s\n", comment_syntax ? comment_syntax : "NULL");

    if (comment_syntax) {
        gtk_text_buffer_insert(buffer, &iter, comment_syntax, -1); // Insert at cursor position
        g_print("Comment inserted at cursor position successfully\n");
    } else {
        g_print("No matching format for comment insertion\n");
    }

    gtk_text_view_place_cursor_onscreen(GTK_TEXT_VIEW(text_view));
    gtk_widget_grab_focus(GTK_WIDGET(text_view));
}
//END COMMENT

///PASTE IMAGE
// Function to decode Base64
static guchar *base64_decode(const gchar *data, gsize *length) {
    return g_base64_decode(data, length);
}

// Function to insert an image into a text view
static gboolean paste_image(GtkWidget *text_view) {
    GtkClipboard *clipboard = gtk_widget_get_clipboard(text_view, GDK_SELECTION_CLIPBOARD);
    GdkPixbuf *pixbuf = gtk_clipboard_wait_for_image(clipboard);

    if (pixbuf) {
        GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
        GtkTextIter iter;
        gtk_text_buffer_get_end_iter(buffer, &iter);
        gtk_text_buffer_insert_pixbuf(buffer, &iter, pixbuf);
        g_object_unref(pixbuf);

        // Clear the clipboard to avoid pasting the link
        gtk_clipboard_clear(clipboard);
        gtk_clipboard_set_text(clipboard, "", -1);

        return TRUE; // Handle the image paste and suppress default paste behavior
    }

    gchar *text = gtk_clipboard_wait_for_text(clipboard);

    if (text && g_str_has_prefix(text, "data:image")) {
        // Handle Base64 image data
        gchar *comma = strchr(text, ',');
        if (comma) {
            gsize length;
            guchar *image_data = base64_decode(comma + 1, &length);
            GInputStream *stream = g_memory_input_stream_new_from_data(image_data, length, g_free);
            GdkPixbuf *pixbuf = gdk_pixbuf_new_from_stream(stream, NULL, NULL);

            if (pixbuf) {
                GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
                GtkTextIter iter;
                gtk_text_buffer_get_end_iter(buffer, &iter);
                gtk_text_buffer_insert_pixbuf(buffer, &iter, pixbuf);
                g_object_unref(pixbuf);
            }

            g_object_unref(stream);
            g_free(text);
            gtk_clipboard_clear(clipboard);
            gtk_clipboard_set_text(clipboard, "", -1); // Ensure clipboard is cleared

            return TRUE; // Handle the image paste and suppress default paste behavior
        }
    }

    g_free(text);
    return FALSE; // If not handled, allow default paste behavior
}

// Function to handle paste-clipboard event
static gboolean on_paste_clipboard(GtkWidget *text_view, GdkEvent *event, gpointer data) {
    if (paste_image(text_view)) {
        return TRUE; // Suppress default paste behavior if image was handled
    } else {
        GtkClipboard *clipboard = gtk_widget_get_clipboard(text_view, GDK_SELECTION_CLIPBOARD);
        gtk_text_buffer_paste_clipboard(gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view)), clipboard, NULL, TRUE);
        return TRUE;
    }
}
//END PAST IMAGES
/////////////////END NOTEBOOK

///DOWNLOAD
// Progress update callback to update the progress bar
static void update_progress_bar(WebKitDownload *download, GParamSpec *pspec, gpointer user_data) {
    GtkProgressBar *progress_bar = GTK_PROGRESS_BAR(user_data);
    gdouble progress = webkit_download_get_estimated_progress(download);

 
    // Update the progress bar fraction
    gtk_progress_bar_set_fraction(progress_bar, progress);
    
    // Optional: Show progress as a percentage on the progress bar
    gchar *progress_text = g_strdup_printf("%.2f%%", progress * 100.0);
    gtk_progress_bar_set_text(progress_bar, progress_text);
    g_free(progress_text);

    //g_print("Download Progress: %.2f%%\n", progress * 100.0);
}
void download_started_callback(WebKitDownload *download, gpointer user_data) {
    // Validate and cast the download object
    download = WEBKIT_DOWNLOAD(user_data);
    if (!WEBKIT_IS_DOWNLOAD(download)) {
        g_warning("Failed to cast the received object to WebKitDownload!");
        return;
    }

    g_print("Download started! Download instance: %p\n", (void *)download);

    // Create a GtkDialog
    GtkWidget *download_dialog = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(download_dialog), "Download"); // Set dialog title
    gtk_widget_set_name(download_dialog, "download_dialog"); // Assign a CSS class to the dialog

    // Get the content area of the dialog
    GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(download_dialog));

    // Create a vertical box to organize widgets
    GtkWidget *download_event_entry_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
    gtk_widget_set_name(download_event_entry_box, "download_event_entry"); // Assign CSS class
    gtk_box_pack_start(GTK_BOX(content_area), download_event_entry_box, TRUE, TRUE, 0);

    // Add a progress bar to the dialog
    GtkWidget *progress_bar = gtk_progress_bar_new();
    gtk_widget_set_name(progress_bar, "progress_bar"); 
    gtk_box_pack_start(GTK_BOX(download_event_entry_box), progress_bar, FALSE, FALSE, 0);

    // Add the filename label (initially showing the filename as text)
    GtkWidget *filename_label = gtk_label_new(NULL);
    gtk_widget_set_name(filename_label, "filename_label");
    gtk_box_pack_start(GTK_BOX(download_event_entry_box), filename_label, FALSE, FALSE, 0);

    // Retrieve the filename from multiple fallback methods
    const gchar *filename = NULL;

    // Attempt to get the filename from the response
    WebKitURIResponse *response = webkit_download_get_response(download);
    filename = response ? webkit_uri_response_get_suggested_filename(response) : NULL;

    // Fallback: Extract the filename from the download URI if no filename is suggested
    if (!filename) {
        const gchar *download_uri = webkit_uri_request_get_uri(webkit_download_get_request(download));
        if (download_uri) {
            GFile *file = g_file_new_for_uri(download_uri);
            filename = g_file_get_basename(file);
            g_object_unref(file);
        } else {
            filename = "Unknown Filename"; // Final fallback
        }
    }

    g_print("Filename: %s\n", filename);

    // Update the filename label
    gchar *filename_text = g_strdup_printf("Filename: %s", filename);
    gtk_label_set_text(GTK_LABEL(filename_label), filename_text);
    g_free(filename_text);

    // Show all widgets in the dialog
    gtk_widget_show_all(download_dialog);
    

    // Connect the progress_updated callback to update the progress bar
    g_signal_connect(download, "notify::estimated-progress", G_CALLBACK(update_progress_bar), progress_bar);

    // Connect the "finished" signal to update the filename label with download status
    g_signal_connect(download, "finished", G_CALLBACK(load_finished), filename_label);

    g_print("Connected 'finished' signal and progress monitoring.\n");
}

static void load_finished(WebKitDownload *download, gpointer user_data) {
    GtkWidget *filename_label = GTK_WIDGET(user_data);

    // Attempt to get the filename from the response
    const gchar *filename = NULL;
    WebKitURIResponse *response = webkit_download_get_response(download);
    filename = response ? webkit_uri_response_get_suggested_filename(response) : NULL;

    // Fallback: Extract the filename from the download URI if no filename is suggested
    if (!filename) {
        const gchar *download_uri = webkit_uri_request_get_uri(webkit_download_get_request(download));
        if (download_uri) {
            GFile *file = g_file_new_for_uri(download_uri);
            filename = g_file_get_basename(file);
            g_object_unref(file);
        } else {
            filename = "Unknown Filename"; // Final fallback
        }
    }

    // Get total data received
    guint64 received_data_length = webkit_download_get_received_data_length(download);

    // Check for errors (no data received)
    if (received_data_length == 0) {
        gchar *error_text = g_strdup_printf("Filename: %s\nDownload error!\nNo data received.", filename);
        gtk_label_set_text(GTK_LABEL(filename_label), error_text);
        g_free(error_text);
        g_warning("Download failed! Filename: %s\nNo data received.\n", filename);
        return; // Exit early since it's an error
    }

    // Convert bytes to megabytes
    double received_data_mb = (double)received_data_length / 1048576.0;

    // Update the filename label with the final status and total data received
    gchar *updated_text = g_strdup_printf("Filename: %s\nDownload finished!\nTotal Data Received: %.2f MB", 
                                           filename, received_data_mb);
    gtk_label_set_text(GTK_LABEL(filename_label), updated_text);
    g_free(updated_text);

    g_print("Download finished!\nFilename: %s\nTotal Data Received: %.2f MB\n", filename, received_data_mb);
}

////
void open_download_folder(GtkWidget *widget, GdkEventButton *event, gpointer data) {
    gchar *download_dir = NULL;
    gchar *command = NULL;

#ifdef _WIN32
    // Get the Downloads folder from Windows environment variable
    download_dir = g_strdup_printf("%s\\Downloads", g_getenv("USERPROFILE"));
    command = g_strconcat("explorer \"", download_dir, "\"", NULL);
#else
    // Get the correct Downloads folder path using xdg-user-dir
    download_dir = g_strdup(g_getenv("XDG_DOWNLOAD_DIR"));

    if (!download_dir) {
        GError *error = NULL;
        gchar *command_output = NULL;
        g_spawn_command_line_sync("xdg-user-dir DOWNLOAD", &command_output, NULL, NULL, &error);

        if (!error) {
            download_dir = g_strdup(g_strstrip(command_output));
        }

        g_free(command_output);
    }

    if (download_dir) {
        command = g_strconcat("xdg-open ", download_dir, NULL);
    }
#endif

    if (command) {
        g_spawn_command_line_async(command, NULL);
        g_free(command);
    }

    g_free(download_dir);
}

void open_file(GtkWidget *widget, GdkEventButton *event, gpointer data) {
    const gchar *file_path = (const gchar *)data;

    // Construct the command to open the file
    gchar *command = g_strconcat("xdg-open ", file_path, NULL);

    // Execute the command asynchronously
    g_spawn_command_line_async(command, NULL);

    g_free(command);
}

//DOWNLOAD LIST
GtkWidget *file_list;

// Function to display files in the Downloads folder


// Function to get the file size
gchar* get_file_size(const gchar *file_path) {
    struct stat st;
    if (stat(file_path, &st) == 0) {
        gchar *size = g_format_size(st.st_size);
        return size;
    } else {
        return g_strdup("Unknown");
    }
}

// Function to get the file format
gchar* get_file_format(const gchar *file_name) {
    gchar **split_name = g_strsplit(file_name, ".", -1);
    gchar *format = g_strdup(split_name[g_strv_length(split_name) - 1]);
    g_strfreev(split_name);
    return format;
}

// Function to get the file modification time
gchar* get_file_mod_time(const gchar *file_path) {
    struct stat st;
    if (stat(file_path, &st) == 0) {
        gchar mod_time_str[20];
        strftime(mod_time_str, sizeof(mod_time_str), "%Y-%m-%d %H:%M:%S", localtime(&st.st_mtime));
        return g_strdup(mod_time_str);
    } else {
        return g_strdup("Unknown");
    }
}

//SHOW DOWNLOAD LIST
// Struct to hold file information for sorting
typedef struct {
    gchar *name;
    gchar *path;
    gchar *size;
    gchar *format;
    gchar *mod_time;
    time_t mtime; // For sorting by modification time
} FileEntry;

static gint compare_by_date(gconstpointer a, gconstpointer b) {
    const FileEntry *file_a = (const FileEntry *)a;
    const FileEntry *file_b = (const FileEntry *)b;

    // Sort in descending order (most recent first)
    return file_b->mtime - file_a->mtime;
}

void show_downloads(GtkWidget *file_list) {
    // Clear existing entries
    GList *children = gtk_container_get_children(GTK_CONTAINER(file_list));
    for (GList *iter = children; iter != NULL; iter = g_list_next(iter)) {
        gtk_widget_destroy(GTK_WIDGET(iter->data));
    }
    g_list_free(children);

    gchar *download_dir = NULL;

#ifdef _WIN32
    // Get Downloads folder path on Windows
    download_dir = g_strdup_printf("%s\\Downloads", g_getenv("USERPROFILE"));
#else
    // Get Downloads folder path dynamically on Linux
    download_dir = g_strdup(g_getenv("XDG_DOWNLOAD_DIR"));

    if (!download_dir) {
        GError *error = NULL;
        gchar *command_output = NULL;
        g_spawn_command_line_sync("xdg-user-dir DOWNLOAD", &command_output, NULL, NULL, &error);

        if (!error) {
            download_dir = g_strdup(g_strstrip(command_output));
        }

        g_free(command_output);
    }
#endif

    if (!download_dir) {
        g_print("Failed to get Downloads folder path.\n");
        return;
    }

    // Open the Downloads directory
    DIR *dir = opendir(download_dir);
    if (!dir) {
        g_print("Failed to open Downloads directory: %s\n", download_dir);
        g_free(download_dir);
        return;
    }

    GList *file_entries = NULL;
    struct dirent *entry;
    
    while ((entry = readdir(dir)) != NULL) {
        if (g_strcmp0(entry->d_name, ".") == 0 || g_strcmp0(entry->d_name, "..") == 0) {
            continue;
        }

        gchar *file_path = g_build_filename(download_dir, entry->d_name, NULL);
        struct stat file_stat;

        if (stat(file_path, &file_stat) == -1) {
            g_free(file_path);
            continue;
        }

        // Populate file information
        FileEntry *file_entry = g_new0(FileEntry, 1);
        file_entry->name = g_strdup(entry->d_name);
        file_entry->path = file_path;
        file_entry->size = get_file_size(file_path);
        file_entry->format = get_file_format(entry->d_name);
        file_entry->mod_time = get_file_mod_time(file_path);
        file_entry->mtime = file_stat.st_mtime;

        file_entries = g_list_prepend(file_entries, file_entry);
    }
    closedir(dir);

    // Sort and display files
    file_entries = g_list_sort(file_entries, compare_by_date);

    // Ensure UI elements are properly created and displayed
    for (GList *iter = file_entries; iter != NULL; iter = g_list_next(iter)) {
        FileEntry *file_entry = (FileEntry *)iter->data;

        gchar *file_info = g_strdup_printf("NAME: %s\nSIZE: %s\nFORMAT: %s\nDOWNLOADED: %s", 
                                           file_entry->name, file_entry->size, file_entry->format, file_entry->mod_time);
        GtkWidget *button = gtk_button_new_with_label(file_info);

        GtkWidget *label = gtk_bin_get_child(GTK_BIN(button));
        gtk_widget_set_halign(label, GTK_ALIGN_START);
        gtk_label_set_xalign(GTK_LABEL(label), 0.0);

        gtk_widget_set_name(button, "file_label");

        if (GTK_IS_CONTAINER(file_list)) {
            GtkWidget *row = gtk_list_box_row_new();
            gtk_container_add(GTK_CONTAINER(row), button);
            gtk_list_box_insert(GTK_LIST_BOX(file_list), row, -1);

            g_signal_connect_data(button, "clicked", G_CALLBACK(open_file), file_entry->path, (GClosureNotify)g_free, G_CONNECT_AFTER);
            gtk_widget_set_name(row, "download_row");

            gtk_widget_show_all(row);
        }

        g_free(file_info);
        g_free(file_entry->name);
        g_free(file_entry->size);
        g_free(file_entry->format);
        g_free(file_entry->mod_time);
        g_free(file_entry);
    }

    g_list_free(file_entries);
    g_free(download_dir);
}




// Function to handle download button click
void on_download_button_clicked(GtkWidget *button, GtkWidget *downloadbar) {
    if (gtk_widget_get_visible(downloadbar)) {
        gtk_widget_hide(downloadbar);
    } else {
        gtk_widget_show(downloadbar);
        show_downloads(file_list); // Reload the contents of the download bar
    }
}
///END DOWNLOAD


///COLORPICKER
void on_colorpicker_button_clicked(GtkWidget *button, GtkWidget *colorpickerbar) {
    if (gtk_widget_get_visible(colorpickerbar)) {
        gtk_widget_hide(colorpickerbar);
    } else {
        gtk_widget_show(colorpickerbar);
    }
}

static double brightness = 1.0;

static gboolean on_draw(GtkWidget *widget, cairo_t *cr, gpointer user_data) {
    int width = gtk_widget_get_allocated_width(widget);
    int height = gtk_widget_get_allocated_height(widget);
    int radius = MIN(width, height) / 2;
    int center_x = width / 2;
    int center_y = height / 2;

    // Create a circular mask
    cairo_arc(cr, center_x, center_y, radius, 0, 2 * G_PI);
    cairo_clip(cr);

    // Draw radial gradient
    for (int i = 0; i < 360; i++) {
        double angle = (double)i * G_PI / 180;
        GdkRGBA color;

        // Calculate the base color at each angle
        color.red = cos(angle) / 2 + 0.5;
        color.green = cos(angle - 2 * G_PI / 3) / 2 + 0.5;
        color.blue = cos(angle - 4 * G_PI / 3) / 2 + 0.5;

        // Blend color with black or white based on brightness
        if (brightness >= 1.0) {
            color.red = color.red + (1.0 - color.red) * (brightness - 1.0);
            color.green = color.green + (1.0 - color.green) * (brightness - 1.0);
            color.blue = color.blue + (1.0 - color.blue) * (brightness - 1.0);
        } else {
            color.red = color.red * brightness;
            color.green = color.green * brightness;
            color.blue = color.blue * brightness;
        }

        // Gradient from center to border
        cairo_set_source_rgba(cr, color.red, color.green, color.blue, 1.0);
        cairo_move_to(cr, center_x, center_y);
        cairo_line_to(cr, center_x + radius * cos(angle), center_y + radius * sin(angle));
        cairo_stroke(cr);
    }

// Draw dark gray center
cairo_set_source_rgba(cr, 38.0 / 255.0, 38.0 / 255.0, 38.0 / 255.0, 1.0);
cairo_arc(cr, center_x, center_y, radius / 8, 0, 2 * G_PI); // Adjust radius to make it bigger
cairo_fill(cr);


    // Draw black border
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
    cairo_arc(cr, center_x, center_y, radius, 0, 2 * G_PI);
    cairo_stroke(cr);

    return FALSE;
}

static void on_brightness_changed(GtkRange *range, gpointer user_data) {
    brightness = gtk_range_get_value(range);
    gtk_widget_queue_draw(GTK_WIDGET(user_data));
}

static GdkRGBA selected_color;

// Callback function to draw the selected color in the drawing area
static gboolean drawcolor_callback(GtkWidget *widget, cairo_t *cr, gpointer data) {
    gdk_cairo_set_source_rgba(cr, &selected_color);
    cairo_paint(cr);
    return FALSE;
}

static void on_color_selected(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
    GtkEntry *rgb_entry = GTK_ENTRY(g_object_get_data(G_OBJECT(widget), "rgb_entry"));
    GtkEntry *hex_entry = GTK_ENTRY(g_object_get_data(G_OBJECT(widget), "hex_entry"));
    GtkWidget *showcolor = GTK_WIDGET(g_object_get_data(G_OBJECT(widget), "showcolor"));

    int x = (int)event->x;
    int y = (int)event->y;
    int width = gtk_widget_get_allocated_width(widget);
    int height = gtk_widget_get_allocated_height(widget);
    double center_x = width / 2;
    double center_y = height / 2;

    double dx = x - center_x;
    double dy = y - center_y;
    double angle = atan2(dy, dx);

    selected_color.red = cos(angle) / 2 + 0.5;
    selected_color.green = cos(angle - 2 * G_PI / 3) / 2 + 0.5;
    selected_color.blue = cos(angle - 4 * G_PI / 3) / 2 + 0.5;

    // Apply brightness scaling
    if (brightness >= 1.0) {
        selected_color.red = selected_color.red + (1.0 - selected_color.red) * (brightness - 1.0);
        selected_color.green = selected_color.green + (1.0 - selected_color.green) * (brightness - 1.0);
        selected_color.blue = selected_color.blue + (1.0 - selected_color.blue) * (brightness - 1.0);
    } else {
        selected_color.red = selected_color.red * brightness;
        selected_color.green = selected_color.green * brightness;
        selected_color.blue = selected_color.blue * brightness;
    }

    selected_color.alpha = 1.0;

    char rgb_text[50];
    snprintf(rgb_text, sizeof(rgb_text), "RGB: %.0f, %.0f, %.0f", selected_color.red * 255, selected_color.green * 255, selected_color.blue * 255);
    gtk_entry_set_text(rgb_entry, rgb_text);

    char hex_text[20];
    snprintf(hex_text, sizeof(hex_text), "Hex: #%02X%02X%02X", (int)(selected_color.red * 255), (int)(selected_color.green * 255), (int)(selected_color.blue * 255));
    gtk_entry_set_text(hex_entry, hex_text);

    // Update the drawing area to show the selected color
    gtk_widget_queue_draw(showcolor);
}

static gboolean on_drawbw(GtkWidget *widget, cairo_t *cr, gpointer user_data) {
    int width = gtk_widget_get_allocated_width(widget);
    int height = gtk_widget_get_allocated_height(widget);

    // Draw grayscale gradient
    for (int i = 0; i < width; i++) {
        double ratio = (double)i / width;
        GdkRGBA color;

        color.red = ratio;
        color.green = ratio;
        color.blue = ratio;

        cairo_set_source_rgba(cr, color.red, color.green, color.blue, 1.0);
        cairo_rectangle(cr, i, 0, 1, height);
        cairo_fill(cr);
    }

    return FALSE;
}

// Function to handle color selection and update the drawing area
static void on_colorbw_selected(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
    GtkEntry *rgb_entry = GTK_ENTRY(g_object_get_data(G_OBJECT(widget), "rgb_entry"));
    GtkEntry *hex_entry = GTK_ENTRY(g_object_get_data(G_OBJECT(widget), "hex_entry"));
    GtkWidget *showcolor = GTK_WIDGET(g_object_get_data(G_OBJECT(widget), "showcolor"));
    int x = (int)event->x;
    int width = gtk_widget_get_allocated_width(widget);

    double ratio = (double)x / width;

    selected_color.red = ratio;
    selected_color.green = ratio;
    selected_color.blue = ratio;
    selected_color.alpha = 1.0;

    char rgb_text[50];
    snprintf(rgb_text, sizeof(rgb_text), "RGB: %.0f, %.0f, %.0f", selected_color.red * 255, selected_color.green * 255, selected_color.blue * 255);
    gtk_entry_set_text(rgb_entry, rgb_text);

    char hex_text[20];
    snprintf(hex_text, sizeof(hex_text), "Hex: #%02X%02X%02X", (int)(selected_color.red * 255), (int)(selected_color.green * 255), (int)(selected_color.blue * 255));
    gtk_entry_set_text(hex_entry, hex_text);

    // Update the drawing area to show the selected color
    gtk_widget_queue_draw(showcolor);
    // Debugging statement
    // //g_print("Color selected: red=%.2f, green=%.2f, blue=%.2f, alpha=%.2f\n", selected_color.red, selected_color.green, selected_color.blue, selected_color.alpha);
}

//////////END COLOR PICKER

///CONSOLE
static void console_button_clicked(GtkWidget *widget, gpointer data) {
    g_return_if_fail(GTK_IS_WIDGET(widget)); // Ensure the widget is valid

    GtkNotebook *notebook = GTK_NOTEBOOK(data);
    if (!GTK_IS_NOTEBOOK(notebook)) {
        g_critical("Invalid GtkNotebook pointer. Ensure 'data' is a GtkNotebook.");
        return;
    }

    gint current_page = gtk_notebook_get_current_page(notebook);
    g_message("Current page index: %d", current_page);

    GtkWidget *page_widget = gtk_notebook_get_nth_page(notebook, current_page);
    if (!page_widget) {
        g_critical("Page widget is NULL. Ensure a valid page exists at index %d.", current_page);
        return;
    }

    WebKitWebView *web_view = WEBKIT_WEB_VIEW(page_widget);
    if (!WEBKIT_IS_WEB_VIEW(web_view)) {
        g_critical("Page widget is not a valid WebKitWebView. Ensure the notebook page contains a WebKitWebView.");
        return;
    }

    WebKitWebInspector *inspector = webkit_web_view_get_inspector(web_view);
    if (!WEBKIT_IS_WEB_INSPECTOR(inspector)) {
        g_critical("Failed to retrieve WebKitWebInspector. Ensure the WebKitWebView is valid.");
        return;
    }

    webkit_web_inspector_show(inspector);
    g_message("WebInspector is now shown.");

    WebKitWebViewBase *inspector_web_view = webkit_web_inspector_get_web_view(inspector);
if (!GTK_IS_WIDGET(GTK_WIDGET(inspector_web_view))) {
    g_critical("Failed to retrieve WebInspector's web view.");
    return;
}

    gtk_widget_set_name(GTK_WIDGET(inspector_web_view), "inspector");
    g_message("Inspector web view name set to 'inspector'.");
}

//END CONSOLE

//SPEED TEST
void on_speedtest_button_clicked(GtkWidget *button, GtkWidget *speedtestbar) {
    if (gtk_widget_get_visible(speedtestbar)) {
        gtk_widget_hide(speedtestbar);
    } else {
        gtk_widget_show(speedtestbar);
    }
}

// Dummy write callback for data
static size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
    return size * nmemb; // Discard data
}

// Helper function to format speed to two decimal places
std::string format_speed(double speed) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << speed;
    return out.str();
}

// Function to measure download speed using Gio/Glib
double measure_download_speed(const std::string& url) {
    GError* error = NULL;
    GFile* file = g_file_new_for_uri(url.c_str());
    GFileInputStream* stream = g_file_read(file, NULL, &error);

    if (error) {
        g_error_free(error);
        return -1; // Failure
    }

    char buffer[8192]; // Larger buffer size
    gsize bytes_read;
    size_t total_bytes = 0; // Track total data
    auto start = std::chrono::high_resolution_clock::now();

    // Warm-up phase
    g_input_stream_read(G_INPUT_STREAM(stream), buffer, sizeof(buffer), NULL, &error);

    // Main read phase
    while ((bytes_read = g_input_stream_read(G_INPUT_STREAM(stream), buffer, sizeof(buffer), NULL, &error)) > 0) {
        total_bytes += bytes_read;
    }

    auto end = std::chrono::high_resolution_clock::now();
    g_object_unref(file);
    g_object_unref(stream);

    if (error) {
        g_error_free(error);
        return -1; // Failure
    }

    std::chrono::duration<double> elapsed = end - start;
    if (elapsed.count() > 0) {
        double speed_mbps = (total_bytes * 8) / (elapsed.count() * 1'000'000); // Convert to Mbps
        return std::round(speed_mbps * 100) / 100; // Round to 2 decimals
    } else {
        return -1; // Prevent division by zero
    }
}


double measure_upload_speed(const std::string& server_ip, int server_port, const std::string& endpoint, const std::string& data) {
    struct addrinfo hints{}, *res;
    hints.ai_family = AF_INET; // IPv4
    hints.ai_socktype = SOCK_STREAM; // TCP

    if (getaddrinfo(server_ip.c_str(), std::to_string(server_port).c_str(), &hints, &res) != 0) {
        std::cerr << "Failed to resolve hostname." << std::endl;
        return -1.0;
    }

    int sock = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
    if (sock < 0) {
        std::cerr << "Error creating socket." << std::endl;
        freeaddrinfo(res);
        return -1.0;
    }

    if (connect(sock, res->ai_addr, res->ai_addrlen) < 0) {
        std::cerr << "Connection to server failed." << std::endl;
        close(sock);
        freeaddrinfo(res);
        return -1.0;
    }

    freeaddrinfo(res);

    // Construct the HTTP POST request
    std::string request = "POST " + endpoint + " HTTP/1.1\r\n";
    request += "Host: " + server_ip + "\r\n";
    request += "Content-Type: application/x-www-form-urlencoded\r\n";
    request += "Content-Length: " + std::to_string(data.size()) + "\r\n\r\n";
    request += data;

    // Measure elapsed time during data transmission
    auto start = std::chrono::high_resolution_clock::now();
    size_t bytes_sent = 0;
    while (bytes_sent < request.size()) {
        int sent = send(sock, request.c_str() + bytes_sent, request.size() - bytes_sent, 0);
        if (sent < 0) {
            std::cerr << "Failed to send data to the server." << std::endl;
            close(sock);
            return -1.0;
        }
        bytes_sent += sent;
    }
    auto end = std::chrono::high_resolution_clock::now();

    // Ignore server response timing
    char buffer[8192] = {0};
    recv(sock, buffer, sizeof(buffer), 0);

    close(sock);

    double elapsed_seconds = std::chrono::duration<double>(end - start).count();
    if (elapsed_seconds == 0) {
        std::cerr << "Elapsed time is zero, cannot calculate speed." << std::endl;
        return -1.0;
    }

    // Calculate upload speed
    double upload_speed_mbps = (bytes_sent * 8.0 / elapsed_seconds) / 1'000'000.0;

    // Log debug information
    std::cout << "Data size (bytes): " << data.size() << "\n";
    std::cout << "Bytes sent: " << bytes_sent << "\n";
    std::cout << "Elapsed time (seconds): " << elapsed_seconds << "\n";
    std::cout << "Calculated upload speed (Mbps): " << upload_speed_mbps << "\n";

    return upload_speed_mbps;
}


// Callback to display speed results in GTK entry fields
void on_test_speed_clicked(GtkWidget* button, gpointer user_data) {
    // Unpack user data (array containing download entry, upload entry, and spinner)
    GtkWidget** widgets = (GtkWidget**)user_data;

    GtkEntry* download_entry = GTK_ENTRY(widgets[0]);
    GtkEntry* upload_entry = GTK_ENTRY(widgets[1]);
    GtkSpinner* spinner = GTK_SPINNER(widgets[2]);

    // Test URL details for download speed testing
    std::string server_ip = "proof.ovh.net";    // For download testing
    int server_port = 80;                       // HTTP port
    std::string endpoint = "/files/100Mb.dat";  // Path for download testing
    std::string upload_server_ip = "httpbin.org";  // Server for upload testing
    std::string upload_endpoint = "/post";      // Path for upload testing
    std::string dummy_data = "Test upload data"; // Data to upload

    // Start the spinner for visual feedback
    gtk_spinner_start(spinner);
    gtk_widget_show(GTK_WIDGET(spinner));

    // Run the speed test in a separate thread to avoid blocking the UI
    std::thread speed_test_thread([=]() {
        // Measure download speed
        double download_speed = measure_download_speed("https://" + server_ip + endpoint);

        // Measure upload speed using a separate server and endpoint
        double upload_speed = measure_upload_speed(upload_server_ip, server_port, upload_endpoint, dummy_data);

        // Package the results and widgets into a structure
        struct SpeedTestResult {
            GtkWidget* download_entry;
            GtkWidget* upload_entry;
            GtkWidget* spinner;
            double download_speed;
            double upload_speed;
        };

        SpeedTestResult* result = new SpeedTestResult{
            GTK_WIDGET(download_entry),
            GTK_WIDGET(upload_entry),
            GTK_WIDGET(spinner),
            download_speed,
            upload_speed
        };

        // Use g_idle_add to safely update the GTK UI from the background thread
        g_idle_add([](gpointer data) -> gboolean {
            // Unpack the result
            SpeedTestResult* result = static_cast<SpeedTestResult*>(data);

            GtkEntry* download_entry = GTK_ENTRY(result->download_entry);
            GtkEntry* upload_entry = GTK_ENTRY(result->upload_entry);
            GtkSpinner* spinner = GTK_SPINNER(result->spinner);

            // Update the UI with results
            if (result->download_speed > 0) {
                gtk_entry_set_text(download_entry, format_speed(result->download_speed).c_str());
            } else {
                gtk_entry_set_text(download_entry, "Download failed");
            }

            if (result->upload_speed > 0) {
                gtk_entry_set_text(upload_entry, format_speed(result->upload_speed).c_str());
            } else {
                gtk_entry_set_text(upload_entry, "Upload failed");
            }

            // Stop and hide the spinner
            gtk_spinner_stop(spinner);
            gtk_widget_hide(GTK_WIDGET(spinner));

            // Free the allocated structure
            delete result;

            return FALSE; // Remove the callback after execution
        }, result);
    });

    // Detach the thread to let it run independently
    speed_test_thread.detach();
}

//SEARCHENGINEMENU (link menu)
void on_open_search_button_clicked(GtkWidget *button, GtkWidget *searchgrid) {
    if (gtk_widget_get_visible(searchgrid)) {
        gtk_widget_hide(searchgrid);
    } else {
        gtk_widget_show(searchgrid);
    }
}

void on_open_news_button_clicked(GtkWidget *button, GtkWidget *news_grid) {
    if (gtk_widget_get_visible(news_grid)) {
        gtk_widget_hide(news_grid);
    } else {
        gtk_widget_show(news_grid);
    }
}

void on_open_video_button_clicked(GtkWidget *button, GtkWidget *video_grid) {
    if (gtk_widget_get_visible(video_grid)) {
        gtk_widget_hide(video_grid);
    } else {
        gtk_widget_show(video_grid);
    }
}

void on_open_education_button_clicked(GtkWidget *button, GtkWidget *education_grid) {
    if (gtk_widget_get_visible(education_grid)) {
        gtk_widget_hide(education_grid);
    } else {
        gtk_widget_show(education_grid);
    }
}
// END SEARCHENGINEMENU (link menu)


//EXTERNAL APPS
//GUIDE
static void on_guidebutton_clicked(GtkWidget *guidebutton, gpointer user_data) {
    // Command to open the guide.cpp app
    const char *command = "/usr/bin/guide"; 
    
    // Run the command asynchronously
    if (g_spawn_command_line_async(command, NULL)) {
        g_print("Guide app launched successfully.\n");
    } else {
        g_warning("Failed to launch the Guide app.\n");
    }
}

//ABOUT
static void on_aboutbutton_clicked(GtkWidget *aboutbutton, gpointer user_data) {
    const char *command = "/usr/bin/about"; 
    
    if (g_spawn_command_line_async(command, NULL)) {
        g_print("About app launched successfully.\n");
    } else {
        g_warning("Failed to launch the About app.\n");
    }
}
//

//CALCULATOR
void on_calculator_button_clicked(GtkWidget *button, GtkWidget *calculatorbar) {
    if (gtk_widget_get_visible(calculatorbar)) {
        gtk_widget_hide(calculatorbar);
    } else {
        gtk_widget_show(calculatorbar);
    }
}

//

static void on_window(GtkWidget *widget, gpointer user_data) {
    load_saved_theme();
    load_saved_engine_choice(&selected_search_engine);
}

/////////INIZIO GUI////////////
int main(int argc, char *argv[]) {
    gtk_init(&argc, &argv);
    gst_init(&argc, &argv);
 
/* 
// Load dark_theme.css when the application starts
load_css("/usr/share/debris/css/dark_theme.css");    
*/

/*switch from root to current user
    if (getuid() == 0) {
        struct passwd *pw = getpwnam(getenv("USER"));  // Get the username from the environment
        if (!pw) {
            fprintf(stderr, "Error: Could not detect logged-in user.\n");
            exit(1);
        }

        printf("App mistakenly running as root. Switching to user: %s\n", pw->pw_name);
        setuid(pw->pw_uid);  // Change process UID to the detected user
    }

    // Confirm that we switched correctly
    printf("Now running as user: %d (%s)\n", getuid(), getpwuid(getuid())->pw_name);*/


// Create a new window

//setenv("WEBKIT_DISABLE_COMPOSITING_MODE", "1", 1);

GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
gtk_window_set_title(GTK_WINDOW(window), "BROWS");
gtk_window_set_default_size(GTK_WINDOW(window), 1440, 900);
gtk_window_set_decorated(GTK_WINDOW(window), false);
gtk_widget_add_events(window, GDK_BUTTON_PRESS_MASK | GDK_POINTER_MOTION_MASK);

g_signal_connect(window, "button-press-event", G_CALLBACK(on_buttonmove_press_event), NULL);
g_signal_connect(window, "destroy", G_CALLBACK(gtk_main_quit), NULL);
g_signal_connect(window, "map", G_CALLBACK(on_window), NULL);
/*g_signal_connect(window, "realize", G_CALLBACK(on_start_web_history_file_create), NULL);*/

if (GTK_IS_WIDGET(window) && GTK_IS_WINDOW(window)) {
    //g_print("Window is a valid top-level window\n");
} else {
    //g_print("Error: Window is not a valid top-level window\n");
}
    
// Create a vertical box to hold the address bar, notebook, and button
GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
gtk_widget_set_name(vbox, "vbox");
gtk_container_add(GTK_CONTAINER(window), vbox);


// Create a NEW box to hold the address bar and main box
GtkWidget *searchbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_widget_set_name(searchbox, "searchbox");
gtk_container_add(GTK_CONTAINER(vbox), searchbox);

// HEADBAR
GtkWidget *headbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
gtk_widget_set_size_request(headbar, 1, 40); // Set height to 40 pixels   
gtk_widget_set_name(headbar, "headbar");
gtk_box_pack_start(GTK_BOX(searchbox), headbar, FALSE, TRUE, 0); // Use gtk_box_pack_start for GTK3

// Create a label and pack it in the center
GtkWidget *label = gtk_label_new("debris");
gtk_widget_set_name(label, "debris_label");
gtk_box_pack_start(GTK_BOX(headbar), label, TRUE, TRUE, 0); // Add label to the headbar

// Button container
GtkWidget *button_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
gtk_widget_set_name(button_box, "button_box");
gtk_box_pack_end(GTK_BOX(headbar), button_box, FALSE, TRUE, 0); // Add button box to the headbar

// Close button
GtkWidget *close_button = gtk_button_new_with_label("");
gtk_widget_set_name(close_button, "close_button");
gtk_widget_set_tooltip_text(close_button, "Close");
gtk_widget_set_size_request(close_button, 24, 24);  
g_signal_connect(close_button, "clicked", G_CALLBACK(on_close_button_clicked), window);
gtk_box_pack_end(GTK_BOX(button_box), close_button, FALSE, TRUE, 0);

// Maximize button
GtkWidget *maximize_button = gtk_button_new_with_label("");
gtk_widget_set_name(maximize_button, "maximize_button");
gtk_widget_set_tooltip_text(maximize_button, "Maximize");
gtk_widget_set_size_request(maximize_button, 24, 24);
g_signal_connect(maximize_button, "clicked", G_CALLBACK(on_maximize_button_clicked), window);
gtk_box_pack_end(GTK_BOX(button_box), maximize_button, FALSE, TRUE, 0);

// Minimize button
GtkWidget *minimize_button = gtk_button_new_with_label("");
gtk_widget_set_name(minimize_button, "minimize_button");
gtk_widget_set_tooltip_text(minimize_button, "Minimize");
gtk_widget_set_size_request(minimize_button, 24, 24);
g_signal_connect(minimize_button, "clicked", G_CALLBACK(on_minimize_button_clicked), window);
gtk_box_pack_end(GTK_BOX(button_box), minimize_button, FALSE, TRUE, 0);



GtkWidget *entry = gtk_entry_new();
gtk_widget_set_name(entry, "entry");
gtk_entry_set_text(GTK_ENTRY(entry), "https://"); 
//gtk_widget_hide(entry); // Hide the address bar initially
gtk_box_pack_start(GTK_BOX(searchbox), entry, FALSE, FALSE, 0);

// ADDRESS BAR


// Create a horizontal box to hold the WebView and the toolbar
GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
gtk_widget_set_name(hbox, "hbox");
gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

// Create the GtkNotebook and set it up
GtkNotebook *notebook = GTK_NOTEBOOK(gtk_notebook_new());
gtk_widget_set_name(GTK_WIDGET(notebook), "notebook"); 
gtk_notebook_set_scrollable(notebook, TRUE); 
gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(notebook), TRUE, TRUE, 0);



//MAIN PAGE
GtkWidget *mainpage_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_widget_set_name(mainpage_box, "mainpage_box");
gtk_widget_set_halign(mainpage_box, GTK_ALIGN_FILL);
gtk_widget_set_valign(mainpage_box, GTK_ALIGN_FILL);
//gtk_box_pack_start(GTK_BOX(notebook), mainpage_box, TRUE, TRUE, 0);
//Sgtk_notebook_append_page(notebook, mainpage_box, gtk_label_new("debris"));

// Create and add a top spacer
GtkWidget *top_spacer = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_box_pack_start(GTK_BOX(mainpage_box), top_spacer, TRUE, TRUE, 0);

// Create and add the mainpage_icon button
GtkWidget *mainpage_icon = gtk_button_new_with_label("");
gtk_widget_set_size_request(mainpage_icon, 250, 250);
gtk_widget_set_name(mainpage_icon, "mainpage_icon");
gtk_widget_set_halign(mainpage_icon, GTK_ALIGN_CENTER);
gtk_widget_set_valign(mainpage_icon, GTK_ALIGN_CENTER);
gtk_box_pack_start(GTK_BOX(mainpage_box), mainpage_icon, FALSE, FALSE, 0);

// Create and add the mainpage_grid
GtkWidget *mainpage_grid = gtk_grid_new();
gtk_widget_set_name(mainpage_grid, "mainpage_grid");
gtk_widget_set_halign(mainpage_grid, GTK_ALIGN_CENTER);
gtk_widget_set_valign(mainpage_grid, GTK_ALIGN_CENTER);
gtk_box_pack_start(GTK_BOX(mainpage_box), mainpage_grid, FALSE, FALSE, 0);

GtkWidget *mainpage_entry = gtk_entry_new();
gtk_widget_set_size_request(mainpage_entry, 800, 50);
gtk_widget_set_name(mainpage_entry, "mainpage_entry");
gtk_entry_set_text(GTK_ENTRY(mainpage_entry), "");
gtk_grid_attach(GTK_GRID(mainpage_grid), mainpage_entry, 0, 0, 1, 1);


// Create and add a bottom spacer
GtkWidget *bottom_spacer = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_box_pack_start(GTK_BOX(mainpage_box), bottom_spacer, TRUE, TRUE, 0);
 
    // Create initial tab label
    GtkWidget *initial_tab_label = gtk_label_new("debris");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), mainpage_box, initial_tab_label);
    gtk_widget_show_all(mainpage_box);

    // Connect signals
    g_signal_connect(mainpage_entry, "activate", G_CALLBACK(on_mainpage_entry_activate), notebook);
    g_signal_connect(entry, "activate", G_CALLBACK(on_url_entry_activate), notebook);


// Add a guard condition to prevent errors on first use
if (GTK_IS_ENTRY(entry)) {
   // gtk_entry_set_text(GTK_ENTRY(entry), "https://www.bing.com");
   gtk_entry_set_text(GTK_ENTRY(entry), "https://");
}


//MENUBAR
    // Create a vertical box for the menu
    GtkWidget *menubar = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_widget_set_name(menubar, "menubar");
    gtk_widget_set_size_request(menubar, 200, 1); // Set width to 48 pixels   
    gtk_box_pack_start(GTK_BOX(hbox), menubar, false, true, 0);
    
    // OPEN BUTTON / TEST button
    GtkWidget *openfilebutton = gtk_button_new_with_label("OPEN");
    //gtk_widget_set_name(openfilebutton, "openfilebutton");
    gtk_widget_set_name(openfilebutton, "mainmenu_open_button");
    gtk_widget_set_halign(openfilebutton, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(menubar), openfilebutton, false, false, 0);
    g_signal_connect(openfilebutton, "clicked", G_CALLBACK(on_openfilebutton_clicked), notebook);

    //SEARCH SETTINGS
    GtkWidget *settings_searchbutton = gtk_button_new_with_label("SEARCH");
    gtk_widget_set_name(settings_searchbutton, "mainmenu_button");
    gtk_widget_set_halign(settings_searchbutton, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(menubar), settings_searchbutton, false, false, 0);
    
    GtkWidget *settings_searchbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_widget_set_name(settings_searchbox, "themebox");
    gtk_box_pack_start(GTK_BOX(menubar), settings_searchbox, TRUE, TRUE, 0);

    GtkWidget *settings_search_grid = gtk_grid_new();
    gtk_container_set_border_width(GTK_CONTAINER(settings_search_grid), 0);
    gtk_widget_set_halign(settings_search_grid, GTK_ALIGN_CENTER);
    gtk_box_pack_start(GTK_BOX(settings_searchbox), settings_search_grid, true, true, 0);   
 
    GtkWidget *settings_search_bing = gtk_button_new_with_label("");
    gtk_widget_set_size_request(settings_search_bing, 36, 36);
    gtk_grid_attach(GTK_GRID(settings_search_grid), settings_search_bing, 0, 0, 1, 1);
    gtk_widget_set_name(settings_search_bing, "settings_search_bing");
    gtk_widget_set_tooltip_text(settings_search_bing, "Bing");
    g_signal_connect(settings_search_bing, "clicked", G_CALLBACK(on_settings_search_bing_clicked), notebook);  
    
    GtkWidget *settings_search_duck = gtk_button_new_with_label("");
    gtk_widget_set_size_request(settings_search_duck, 36, 36);
    gtk_grid_attach(GTK_GRID(settings_search_grid), settings_search_duck, 1, 0, 1, 1);
    gtk_widget_set_name(settings_search_duck, "settings_search_duck");
    gtk_widget_set_tooltip_text(settings_search_duck, "DuckDuckGO");
    g_signal_connect(settings_search_duck, "clicked", G_CALLBACK(on_settings_search_duck_clicked), notebook);

    GtkWidget *settings_search_google = gtk_button_new_with_label("");
    gtk_widget_set_size_request(settings_search_google, 36, 36);
    gtk_grid_attach(GTK_GRID(settings_search_grid), settings_search_google, 2, 0, 1, 1);
    gtk_widget_set_name(settings_search_google, "settings_search_google");
    gtk_widget_set_tooltip_text(settings_search_google, "Google");
    g_signal_connect(settings_search_google, "clicked", G_CALLBACK(on_settings_search_google_clicked), notebook);

    GtkWidget *settings_search_yahoo = gtk_button_new_with_label("");
    gtk_widget_set_size_request(settings_search_yahoo, 36, 36);
    gtk_grid_attach(GTK_GRID(settings_search_grid), settings_search_yahoo, 0, 1, 1, 1);
    gtk_widget_set_name(settings_search_yahoo, "settings_search_yahoo");
    gtk_widget_set_tooltip_text(settings_search_yahoo, "Yahoo");
    g_signal_connect(settings_search_yahoo, "clicked", G_CALLBACK(on_settings_search_yahoo_clicked), notebook); 

    GtkWidget *settings_search_wiki = gtk_button_new_with_label("");
    gtk_widget_set_size_request(settings_search_wiki, 36, 36);
    gtk_grid_attach(GTK_GRID(settings_search_grid), settings_search_wiki, 1, 1, 1, 1);
    gtk_widget_set_name(settings_search_wiki, "settings_search_wiki");
    gtk_widget_set_tooltip_text(settings_search_wiki, "Wikipedia");
    g_signal_connect(settings_search_wiki, "clicked", G_CALLBACK(on_settings_search_wiki_clicked), notebook);    
    
    GtkWidget *settings_search_wolframalpha = gtk_button_new_with_label("");
    gtk_widget_set_size_request(settings_search_wolframalpha, 36, 36);
    gtk_grid_attach(GTK_GRID(settings_search_grid), settings_search_wolframalpha, 2, 1, 1, 1);
    gtk_widget_set_name(settings_search_wolframalpha, "settings_search_wolframalpha");
    gtk_widget_set_tooltip_text(settings_search_wolframalpha, "Wolframalpha");
    g_signal_connect(settings_search_wolframalpha, "clicked", G_CALLBACK(on_settings_search_wolframalpha_clicked), notebook);     
    //END SEARCH SETTINGS
    
    //PRIVACY
    GtkWidget *privacybutton = gtk_button_new_with_label("PRIVACY");
    gtk_widget_set_name(privacybutton, "mainmenu_button");
    gtk_widget_set_halign(privacybutton, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(menubar), privacybutton, false, false, 0);

    GtkWidget *privacybox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_widget_set_name(privacybox, "themebox");
    gtk_box_pack_start(GTK_BOX(menubar), privacybox, TRUE, TRUE, 0);

    GtkWidget *privacy_grid = gtk_grid_new();
    gtk_container_set_border_width(GTK_CONTAINER(privacy_grid), 5);
    gtk_widget_set_halign(privacy_grid, GTK_ALIGN_CENTER);
    gtk_box_pack_start(GTK_BOX(privacybox), privacy_grid, TRUE, TRUE, 0);   
  
    GtkWidget *delall_bookmarks_button = gtk_button_new_with_label("");
    gtk_widget_set_name(delall_bookmarks_button, "privacy_button_bookmarks");
    gtk_widget_set_size_request(delall_bookmarks_button, 36, 36);
    gtk_grid_attach(GTK_GRID(privacy_grid), delall_bookmarks_button, 0, 1, 1, 1);
    gtk_widget_set_tooltip_text(delall_bookmarks_button, "Delete all bookmarks");
    g_signal_connect(delall_bookmarks_button, "clicked", G_CALLBACK(on_delall_bookmarks_button_clicked), notebook);

    GtkWidget *delall_history_button = gtk_button_new_with_label("");
    gtk_widget_set_name(delall_history_button, "privacy_button_history");
    gtk_widget_set_size_request(delall_history_button, 36, 36);
    gtk_grid_attach(GTK_GRID(privacy_grid), delall_history_button, 1, 1, 1, 1);
    gtk_widget_set_tooltip_text(delall_history_button, "Delete all history");
    g_signal_connect(delall_history_button, "clicked", G_CALLBACK(on_delall_history_button_clicked), notebook);

    GtkWidget *delall_cookies_button = gtk_button_new_with_label("");
    gtk_widget_set_name(delall_cookies_button, "privacy_button_cookies");
    gtk_widget_set_size_request(delall_cookies_button, 36, 36);
    gtk_grid_attach(GTK_GRID(privacy_grid), delall_cookies_button, 0, 2, 1, 1);
    gtk_widget_set_tooltip_text(delall_cookies_button, "Delete all cookies");
    g_signal_connect(delall_cookies_button, "clicked", G_CALLBACK(on_delall_cookies_button_clicked), notebook);
    
    GtkWidget *delall_caches_button = gtk_button_new_with_label("");
    gtk_widget_set_name(delall_caches_button, "privacy_button_caches");
    gtk_widget_set_size_request(delall_caches_button, 36, 36);
    gtk_grid_attach(GTK_GRID(privacy_grid), delall_caches_button, 1, 2, 1, 1);
    gtk_widget_set_tooltip_text(delall_caches_button, "Delete all caches");
    g_signal_connect(delall_caches_button, "clicked", G_CALLBACK(on_delall_caches_button_clicked), notebook);    
    
    //g_signal_connect(privacybutton, "clicked", G_CALLBACK(on_privacybutton_clicked), privacybox);   
    //END PRIVACY

    //THEME
    GtkWidget *themebutton = gtk_button_new_with_label("THEME");
    gtk_widget_set_name(themebutton, "mainmenu_button");
    gtk_widget_set_halign(themebutton, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(menubar), themebutton, false, false, 0);
    
    GtkWidget *themebox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_widget_set_name(themebox, "themebox");
    gtk_box_pack_start(GTK_BOX(menubar), themebox, TRUE, TRUE, 0);

    GtkWidget *theme_grid = gtk_grid_new();
    gtk_container_set_border_width(GTK_CONTAINER(theme_grid), 5);
    gtk_widget_set_halign(theme_grid, GTK_ALIGN_CENTER);
    gtk_box_pack_start(GTK_BOX(themebox), theme_grid, TRUE, TRUE, 0);   
  
    GtkWidget *themebutton_light = gtk_button_new_with_label("");
    gtk_widget_set_size_request(themebutton_light, 48, 48);
    gtk_grid_attach(GTK_GRID(theme_grid), themebutton_light, 0, 1, 1, 1);
    gtk_widget_set_name(themebutton_light, "button_themelight");
    gtk_widget_set_tooltip_text(themebutton_light, "Light");
    g_signal_connect(themebutton_light, "clicked", G_CALLBACK(on_themebutton_light_clicked), NULL);
   
    GtkWidget *themebutton_dark = gtk_button_new_with_label("");
    gtk_widget_set_size_request(themebutton_dark, 48, 48);
    gtk_grid_attach(GTK_GRID(theme_grid), themebutton_dark, 1, 1, 1, 1);
    gtk_widget_set_name(themebutton_dark, "button_themedark");
    gtk_widget_set_tooltip_text(themebutton_dark, "Dark");
    g_signal_connect(themebutton_dark, "clicked", G_CALLBACK(on_themebutton_dark_clicked), NULL);
    
    GtkWidget *themebutton_darkcolor = gtk_button_new_with_label("");
    gtk_widget_set_size_request(themebutton_darkcolor, 48, 48);
    gtk_grid_attach(GTK_GRID(theme_grid), themebutton_darkcolor, 0, 2, 1, 1);
    gtk_widget_set_name(themebutton_darkcolor, "button_themedarkbow");
    gtk_widget_set_tooltip_text(themebutton_darkcolor, "DarkBow");
    g_signal_connect(themebutton_darkcolor, "clicked", G_CALLBACK(on_themebutton_darkcolor_clicked), NULL);

    GtkWidget *themebutton_lightcolor = gtk_button_new_with_label("");
    gtk_widget_set_size_request(themebutton_lightcolor, 48, 48);
    gtk_grid_attach(GTK_GRID(theme_grid), themebutton_lightcolor, 1, 2, 1, 1);
    gtk_widget_set_name(themebutton_lightcolor, "button_themelightbow");
    gtk_widget_set_tooltip_text(themebutton_lightcolor, "LightBow");
    g_signal_connect(themebutton_lightcolor, "clicked", G_CALLBACK(on_themebutton_lightcolor_clicked), NULL);

    GtkWidget *themebutton_matrix = gtk_button_new_with_label("");
    gtk_widget_set_size_request(themebutton_matrix, 48, 48);
    gtk_grid_attach(GTK_GRID(theme_grid), themebutton_matrix, 0, 3, 1, 1);
    gtk_widget_set_name(themebutton_matrix, "button_themematrix");
    gtk_widget_set_tooltip_text(themebutton_matrix, "Matrix");
    g_signal_connect(themebutton_matrix, "clicked", G_CALLBACK(on_themebutton_matrix_clicked), NULL);

    GtkWidget *insidemenu_accesstheme_title = gtk_button_new_with_label("Accessibility theme");
    gtk_widget_set_name(insidemenu_accesstheme_title, "mainmenu_button");
    gtk_grid_attach(GTK_GRID(theme_grid), insidemenu_accesstheme_title, 0, 4, 2, 1);
    
    GtkWidget *themebutton_accessibility1 = gtk_button_new_with_label("");
    gtk_widget_set_size_request(themebutton_accessibility1, 48, 48);
    gtk_grid_attach(GTK_GRID(theme_grid), themebutton_accessibility1, 0, 5, 1, 1);
    gtk_widget_set_name(themebutton_accessibility1, "button_theme_accessibility1");
    gtk_widget_set_tooltip_text(themebutton_accessibility1, "Accessibility light");
    g_signal_connect(themebutton_accessibility1, "clicked", G_CALLBACK(on_themebutton_accessibility1_clicked), NULL);
    
    GtkWidget *themebutton_accessibility_dark = gtk_button_new_with_label("");
    gtk_widget_set_size_request(themebutton_accessibility_dark, 48, 48);
    gtk_grid_attach(GTK_GRID(theme_grid), themebutton_accessibility_dark, 1, 5, 1, 1);
    gtk_widget_set_name(themebutton_accessibility_dark, "button_theme_accessibility_dark");
    gtk_widget_set_tooltip_text(themebutton_accessibility_dark, "Accessibility dark");
    g_signal_connect(themebutton_accessibility_dark, "clicked", G_CALLBACK(on_themebutton_accessibility_dark_clicked), NULL);

    GtkWidget *themebutton_accessibility_dark_cyan = gtk_button_new_with_label("");
    gtk_widget_set_size_request(themebutton_accessibility_dark_cyan, 48, 48);
    gtk_grid_attach(GTK_GRID(theme_grid), themebutton_accessibility_dark_cyan, 0, 6, 1, 1);
    gtk_widget_set_name(themebutton_accessibility_dark_cyan, "button_theme_accessibility_dark_cyan");
    gtk_widget_set_tooltip_text(themebutton_accessibility_dark_cyan, "Accessibility Black Contrast");
    g_signal_connect(themebutton_accessibility_dark_cyan, "clicked", G_CALLBACK(on_themebutton_accessibility_dark_cyan_clicked), NULL);
    
    //g_signal_connect(themebutton, "clicked", G_CALLBACK(on_themebutton_clicked), themebox);   
//END THEME

//FAQ
    GtkWidget *guidebutton = gtk_button_new_with_label("GUIDE");
    gtk_widget_set_name(guidebutton, "mainmenu_open_button");
    gtk_widget_set_halign(guidebutton, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(menubar), guidebutton, FALSE, TRUE, 0);
    g_signal_connect(guidebutton, "clicked", G_CALLBACK(on_guidebutton_clicked), NULL);    

    GtkWidget *aboutbutton = gtk_button_new_with_label("ABOUT");
    gtk_widget_set_name(aboutbutton, "mainmenu_open_button");
    gtk_widget_set_halign(aboutbutton, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(menubar), aboutbutton, FALSE, TRUE, 0);
    g_signal_connect(aboutbutton, "clicked", G_CALLBACK(on_aboutbutton_clicked), NULL);  
//END FAQ

//SEARCHENGINEBAR links
    GtkWidget *searchenginebar = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_widget_set_name(searchenginebar, "searchenginebar");
    gtk_widget_set_size_request(searchenginebar, 220, 1); 
    gtk_box_pack_start(GTK_BOX(hbox), searchenginebar, false, true, 0);

GtkWidget *scrolled_links = gtk_scrolled_window_new(NULL, NULL);
gtk_widget_set_name(scrolled_links, "scrolled_window");
gtk_box_pack_start(GTK_BOX(searchenginebar), scrolled_links, TRUE, TRUE, 0);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_links), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

    GtkWidget *container_links = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_container_add(GTK_CONTAINER(scrolled_links), container_links); 

    GtkWidget *open_search = gtk_button_new_with_label("SEARCH");
    gtk_widget_set_name(open_search, "mainmenu_open_button");
    gtk_widget_set_halign(open_search, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(container_links), open_search, false, false, 0);
    
    GtkWidget *searchgrid = gtk_grid_new();
    gtk_container_set_border_width(GTK_CONTAINER(searchgrid), 5);
    gtk_widget_set_name(searchgrid, "links_grid");    
    gtk_box_pack_start(GTK_BOX(container_links), searchgrid, false, false, 0);
    g_signal_connect(open_search, "clicked", G_CALLBACK(on_open_search_button_clicked), searchgrid);
    
    // BING button
    GtkWidget *bingbutton = gtk_button_new_with_label("");
    gtk_widget_set_name(bingbutton, "bing_button");
    gtk_widget_set_tooltip_text(bingbutton, "Bing");
    gtk_grid_attach(GTK_GRID(searchgrid), bingbutton, 0, 0, 1, 1);
    g_signal_connect(bingbutton, "clicked", G_CALLBACK(on_bing_button_clicked), notebook);  
    // DUCK button
    GtkWidget *duckbutton = gtk_button_new_with_label("");
    gtk_widget_set_name(duckbutton, "duck_button");
    gtk_widget_set_tooltip_text(duckbutton, "DuckDuckGo");
    gtk_grid_attach(GTK_GRID(searchgrid), duckbutton, 1, 0, 1, 1);
    g_signal_connect(duckbutton, "clicked", G_CALLBACK(on_duck_button_clicked), notebook);  
    // GOOGLE button
    GtkWidget *googlebutton = gtk_button_new_with_label("");
    gtk_widget_set_name(googlebutton, "google_button");
    gtk_widget_set_tooltip_text(googlebutton, "Google");
    gtk_grid_attach(GTK_GRID(searchgrid), googlebutton, 2, 0, 1, 1);
    g_signal_connect(googlebutton, "clicked", G_CALLBACK(on_google_button_clicked), notebook);   
    // WALFRAM button
    GtkWidget *wolframalpha_button = gtk_button_new_with_label("");
    gtk_widget_set_name(wolframalpha_button, "wolfram_button");
    gtk_widget_set_tooltip_text(wolframalpha_button, "WolframAlpha");
    gtk_grid_attach(GTK_GRID(searchgrid), wolframalpha_button, 0, 1, 1, 1);
    g_signal_connect(wolframalpha_button, "clicked", G_CALLBACK(on_wolframalpha_button_clicked), notebook);  
    // WIKI button
    GtkWidget *wikibutton = gtk_button_new_with_label("");
    gtk_widget_set_name(wikibutton, "wiki_button");
    gtk_widget_set_tooltip_text(wikibutton, "Wikipedia");
    gtk_grid_attach(GTK_GRID(searchgrid), wikibutton, 1, 1, 1, 1);
    g_signal_connect(wikibutton, "clicked", G_CALLBACK(on_wiki_button_clicked), notebook);      
    // YAHOO button
    GtkWidget *yahoo_button = gtk_button_new_with_label("");
    gtk_widget_set_name(yahoo_button, "yahoo_button");
    gtk_widget_set_tooltip_text(yahoo_button, "Yahoo");
    gtk_grid_attach(GTK_GRID(searchgrid), yahoo_button, 2, 1, 1, 1);   
    g_signal_connect(yahoo_button, "clicked", G_CALLBACK(on_yahoo_button_clicked), notebook);  
//END SEARCH ENGINE 
//START NEWS
    GtkWidget *open_news = gtk_button_new_with_label("NEWS");
    gtk_widget_set_name(open_news, "mainmenu_open_button");
    gtk_widget_set_halign(open_news, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(container_links), open_news, FALSE, FALSE, 0);
    
    GtkWidget *news_grid = gtk_grid_new();
    gtk_container_set_border_width(GTK_CONTAINER(news_grid), 5);
    gtk_widget_set_name(news_grid, "links_grid");    
    gtk_box_pack_start(GTK_BOX(container_links), news_grid, FALSE, FALSE, 0);
    g_signal_connect(open_news, "clicked", G_CALLBACK(on_open_news_button_clicked), news_grid);
    
    GtkWidget *ap_button = gtk_button_new_with_label("");
    gtk_widget_set_name(ap_button, "ap_button");
    gtk_widget_set_tooltip_text(ap_button, "AP News");
    gtk_grid_attach(GTK_GRID(news_grid), ap_button, 0, 0, 1, 1);   
    g_signal_connect(ap_button, "clicked", G_CALLBACK(on_ap_button_clicked), notebook); 

    GtkWidget *alarabya_button = gtk_button_new_with_label("");
    gtk_widget_set_name(alarabya_button, "alarabya_button");
    gtk_widget_set_tooltip_text(alarabya_button, "Al Arabya");
    gtk_grid_attach(GTK_GRID(news_grid), alarabya_button, 1, 0, 1, 1);   
    g_signal_connect(alarabya_button, "clicked", G_CALLBACK(on_alarabya_button_clicked), notebook); 

    GtkWidget *aljazeera_button = gtk_button_new_with_label("");
    gtk_widget_set_name(aljazeera_button, "aljazeera_button");
    gtk_widget_set_tooltip_text(aljazeera_button, "Aljazeera");
    gtk_grid_attach(GTK_GRID(news_grid), aljazeera_button, 2, 0, 1, 1);   
    g_signal_connect(aljazeera_button, "clicked", G_CALLBACK(on_aljazeera_button_clicked), notebook); 

    GtkWidget *cnbc_button = gtk_button_new_with_label("");
    gtk_widget_set_name(cnbc_button, "cnbc_button");
    gtk_widget_set_tooltip_text(cnbc_button, "CNBC");
    gtk_grid_attach(GTK_GRID(news_grid), cnbc_button, 0, 1, 1, 1);   
    g_signal_connect(cnbc_button, "clicked", G_CALLBACK(on_cnbc_button_clicked), notebook); 

    GtkWidget *lemonde_button = gtk_button_new_with_label("");
    gtk_widget_set_name(lemonde_button, "lemonde_button");
    gtk_widget_set_tooltip_text(lemonde_button, "Le Monde");
    gtk_grid_attach(GTK_GRID(news_grid), lemonde_button, 1, 1, 1, 1);   
    g_signal_connect(lemonde_button, "clicked", G_CALLBACK(on_lemonde_button_clicked), notebook); 

    GtkWidget *ecns_button = gtk_button_new_with_label("");
    gtk_widget_set_name(ecns_button, "ecns_button");
    gtk_widget_set_tooltip_text(ecns_button, "ECNS");
    gtk_grid_attach(GTK_GRID(news_grid), ecns_button, 2, 1, 1, 1);   
    g_signal_connect(ecns_button, "clicked", G_CALLBACK(on_ecns_button_clicked), notebook); 

    GtkWidget *euronews_button = gtk_button_new_with_label("");
    gtk_widget_set_name(euronews_button, "euronews_button");
    gtk_widget_set_tooltip_text(euronews_button, "Euronews");
    gtk_grid_attach(GTK_GRID(news_grid), euronews_button, 0, 2, 1, 1);   
    g_signal_connect(euronews_button, "clicked", G_CALLBACK(on_euronews_button_clicked), notebook); 

    GtkWidget *nature_button = gtk_button_new_with_label("");
    gtk_widget_set_name(nature_button, "nature_button");
    gtk_widget_set_tooltip_text(nature_button, "Nature");
    gtk_grid_attach(GTK_GRID(news_grid), nature_button, 1, 2, 1, 1);   
    g_signal_connect(nature_button, "clicked", G_CALLBACK(on_nature_button_clicked), notebook);  

    GtkWidget *nytimes_button = gtk_button_new_with_label("");
    gtk_widget_set_name(nytimes_button, "nytimes_button");
    gtk_widget_set_tooltip_text(nytimes_button, "New York Times");
    gtk_grid_attach(GTK_GRID(news_grid), nytimes_button, 2, 2, 1, 1);   
    g_signal_connect(nytimes_button, "clicked", G_CALLBACK(on_nytimes_button_clicked), notebook);  

    GtkWidget *moscowtimes_button = gtk_button_new_with_label("");
    gtk_widget_set_name(moscowtimes_button, "moscowtimes_button");
    gtk_widget_set_tooltip_text(moscowtimes_button, "The Moscow Times");
    gtk_grid_attach(GTK_GRID(news_grid), moscowtimes_button, 0, 3, 1, 1);   
    g_signal_connect(moscowtimes_button, "clicked", G_CALLBACK(on_moscowtimes_button_clicked), notebook);  

    GtkWidget *phys_button = gtk_button_new_with_label("");
    gtk_widget_set_name(phys_button, "phys_button");
    gtk_widget_set_tooltip_text(phys_button, "Phys.org");
    gtk_grid_attach(GTK_GRID(news_grid), phys_button, 1, 3, 1, 1);   
    g_signal_connect(phys_button, "clicked", G_CALLBACK(on_phys_button_clicked), notebook);    
    
    GtkWidget *sciencenews_button = gtk_button_new_with_label("");
    gtk_widget_set_name(sciencenews_button, "sciencenews_button");
    gtk_widget_set_tooltip_text(sciencenews_button, "Science News");
    gtk_grid_attach(GTK_GRID(news_grid), sciencenews_button, 2, 3, 1, 1);   
    g_signal_connect(sciencenews_button, "clicked", G_CALLBACK(on_sciencenews_button_clicked), notebook);      
//END NEWS
//START VIDEO
    GtkWidget *open_video = gtk_button_new_with_label("VIDEO");
    gtk_widget_set_name(open_video, "mainmenu_open_button");
    gtk_widget_set_halign(open_video, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(container_links), open_video, FALSE, TRUE, 0);
    
    GtkWidget *video_grid = gtk_grid_new();
    gtk_container_set_border_width(GTK_CONTAINER(video_grid), 5);
    gtk_widget_set_name(video_grid, "links_grid");        
    gtk_box_pack_start(GTK_BOX(container_links), video_grid, FALSE, FALSE, 0);
    g_signal_connect(open_video, "clicked", G_CALLBACK(on_open_video_button_clicked), video_grid);   
    
    // ARTE button
    GtkWidget *arte_button = gtk_button_new_with_label("");
    gtk_widget_set_name(arte_button, "arte_button");
    gtk_widget_set_tooltip_text(arte_button, "Arte");
    gtk_grid_attach(GTK_GRID(video_grid), arte_button, 0, 0, 1, 1);
    g_signal_connect(arte_button, "clicked", G_CALLBACK(on_arte_button_clicked), notebook);  
    // BBC button
    GtkWidget *bbc_button = gtk_button_new_with_label("");
    gtk_widget_set_name(bbc_button, "bbc_button");
    gtk_widget_set_tooltip_text(bbc_button, "BBC");
    gtk_grid_attach(GTK_GRID(video_grid), bbc_button, 1, 0, 1, 1);
    g_signal_connect(bbc_button, "clicked", G_CALLBACK(on_bbc_button_clicked), notebook);  
    // DAiLYMOTION button
    GtkWidget *dailymotion_button = gtk_button_new_with_label("");
    gtk_widget_set_name(dailymotion_button, "dailymotion_button");
    gtk_widget_set_tooltip_text(dailymotion_button, "Dailymotion");
    gtk_grid_attach(GTK_GRID(video_grid), dailymotion_button, 2, 0, 1, 1);
    g_signal_connect(dailymotion_button, "clicked", G_CALLBACK(on_dailymotion_button_clicked), notebook);  
    // DOC+ button
    GtkWidget *docplus_button = gtk_button_new_with_label("");
    gtk_widget_set_name(docplus_button, "docplus_button");
    gtk_widget_set_tooltip_text(docplus_button, "Documentary+");
    gtk_grid_attach(GTK_GRID(video_grid), docplus_button, 0, 1, 1, 1);
    g_signal_connect(docplus_button, "clicked", G_CALLBACK(on_docplus_button_clicked), notebook);  
    // JOVE button
    GtkWidget *jove_button = gtk_button_new_with_label("");
    gtk_widget_set_name(jove_button, "jove_button");
    gtk_widget_set_tooltip_text(jove_button, "Jove+");
    gtk_grid_attach(GTK_GRID(video_grid), jove_button, 1, 1, 1, 1);
    g_signal_connect(jove_button, "clicked", G_CALLBACK(on_jove_button_clicked), notebook);  
    // NASATV button
    GtkWidget *nasatv_button = gtk_button_new_with_label("");
    gtk_widget_set_name(nasatv_button, "nasatv_button");
    gtk_widget_set_tooltip_text(nasatv_button, "Nasa+");
    gtk_grid_attach(GTK_GRID(video_grid), nasatv_button, 2, 1, 1, 1);
    g_signal_connect(nasatv_button, "clicked", G_CALLBACK(on_nasatv_button_clicked), notebook);  
    // ODYSEE button
    GtkWidget *odysee_button = gtk_button_new_with_label("");
    gtk_widget_set_name(odysee_button, "odysee_button");
    gtk_widget_set_tooltip_text(odysee_button, "Odysee");
    gtk_grid_attach(GTK_GRID(video_grid), odysee_button, 0, 2, 1, 1);
    g_signal_connect(odysee_button, "clicked", G_CALLBACK(on_odysee_button_clicked), notebook);  
    // PBS button
    GtkWidget *pbs_button = gtk_button_new_with_label("");
    gtk_widget_set_name(pbs_button, "pbs_button");
    gtk_widget_set_tooltip_text(pbs_button, "PBS");
    gtk_grid_attach(GTK_GRID(video_grid), pbs_button, 1, 2, 1, 1);
    g_signal_connect(pbs_button, "clicked", G_CALLBACK(on_pbs_button_clicked), notebook);  
    // SCIENCE_AAAS button
    GtkWidget *scienceaaas_button = gtk_button_new_with_label("");
    gtk_widget_set_name(scienceaaas_button, "scienceaaas_button");
    gtk_widget_set_tooltip_text(scienceaaas_button, "Science");
    gtk_grid_attach(GTK_GRID(video_grid), scienceaaas_button, 2, 2, 1, 1);
    g_signal_connect(scienceaaas_button, "clicked", G_CALLBACK(on_scienceaaas_button_clicked), notebook);  
    // TED button
    GtkWidget *ted_button = gtk_button_new_with_label("");
    gtk_widget_set_name(ted_button, "ted_button");
    gtk_widget_set_tooltip_text(ted_button, "TED");
    gtk_grid_attach(GTK_GRID(video_grid), ted_button, 0, 3, 1, 1);
    g_signal_connect(ted_button, "clicked", G_CALLBACK(on_ted_button_clicked), notebook);  
    // TWITCH button
    GtkWidget *twitch_button = gtk_button_new_with_label("");
    gtk_widget_set_name(twitch_button, "twitch_button");
    gtk_widget_set_tooltip_text(twitch_button, "Twitch");
    gtk_grid_attach(GTK_GRID(video_grid), twitch_button, 1, 3, 1, 1);
    g_signal_connect(twitch_button, "clicked", G_CALLBACK(on_twitch_button_clicked), notebook);  
    // YOUTUBE button
    GtkWidget *youtube_button = gtk_button_new_with_label("");
    gtk_widget_set_name(youtube_button, "youtube_button");
    gtk_widget_set_tooltip_text(youtube_button, "YouTube");
    gtk_grid_attach(GTK_GRID(video_grid), youtube_button, 2, 3, 1, 1);
    g_signal_connect(youtube_button, "clicked", G_CALLBACK(on_youtube_button_clicked), notebook);  

//END VIDEO
//START EDUCATION
    GtkWidget *open_education = gtk_button_new_with_label("EDUCATION");
    gtk_widget_set_name(open_education, "mainmenu_open_button");
    gtk_widget_set_halign(open_education, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(container_links), open_education, FALSE, TRUE, 0);  
    
    GtkWidget *education_grid = gtk_grid_new();
    gtk_container_set_border_width(GTK_CONTAINER(education_grid), 5);
    gtk_widget_set_name(education_grid, "links_grid");        
    gtk_box_pack_start(GTK_BOX(container_links), education_grid, FALSE, FALSE, 0);
    g_signal_connect(open_education, "clicked", G_CALLBACK(on_open_education_button_clicked), education_grid); 
    
    // BIODIGITAL button
    GtkWidget *biodigital_button = gtk_button_new_with_label("");
    gtk_widget_set_name(biodigital_button, "biodigital_button");
    gtk_widget_set_tooltip_text(biodigital_button, "Biodigital");
    gtk_grid_attach(GTK_GRID(education_grid), biodigital_button, 0, 0, 1, 1);
    g_signal_connect(biodigital_button, "clicked", G_CALLBACK(on_biodigital_button_clicked), notebook);  

    // DLMF button
    GtkWidget *dlmf_button = gtk_button_new_with_label("");
    gtk_widget_set_name(dlmf_button, "dlmf_button");
    gtk_widget_set_tooltip_text(dlmf_button, "Dlmf");
    gtk_grid_attach(GTK_GRID(education_grid), dlmf_button, 1, 0, 1, 1);
    g_signal_connect(dlmf_button, "clicked", G_CALLBACK(on_dlmf_button_clicked), notebook);  
   
    // EFUNDA button
    GtkWidget *efunda_button = gtk_button_new_with_label("");
    gtk_widget_set_name(efunda_button, "efunda_button");
    gtk_widget_set_tooltip_text(efunda_button, "EFunda");
    gtk_grid_attach(GTK_GRID(education_grid), efunda_button, 2, 0, 1, 1);
    g_signal_connect(efunda_button, "clicked", G_CALLBACK(on_efunda_button_clicked), notebook); 
   
    // MIT button
    GtkWidget *mit_button = gtk_button_new_with_label("");
    gtk_widget_set_name(mit_button, "mit_button");
    gtk_widget_set_tooltip_text(mit_button, "MIT");
    gtk_grid_attach(GTK_GRID(education_grid), mit_button, 0, 1, 1, 1);
    g_signal_connect(mit_button, "clicked", G_CALLBACK(on_mit_button_clicked), notebook);    
   
    // ONEZOOM button
    GtkWidget *onezoom_button = gtk_button_new_with_label("");
    gtk_widget_set_name(onezoom_button, "onezoom_button");
    gtk_widget_set_tooltip_text(onezoom_button, "One Zoom");
    gtk_grid_attach(GTK_GRID(education_grid), onezoom_button, 1, 1, 1, 1);
    g_signal_connect(onezoom_button, "clicked", G_CALLBACK(on_onezoom_button_clicked), notebook);     
   
    // PBDB button
    GtkWidget *pbdb_button = gtk_button_new_with_label("");
    gtk_widget_set_name(pbdb_button, "pbdb_button");
    gtk_widget_set_tooltip_text(pbdb_button, "PBDB");
    gtk_grid_attach(GTK_GRID(education_grid), pbdb_button, 2, 1, 1, 1);
    g_signal_connect(pbdb_button, "clicked", G_CALLBACK(on_pbdb_button_clicked), notebook);    
 
    // PT button
    GtkWidget *pt_button = gtk_button_new_with_label("");
    gtk_widget_set_name(pt_button, "pt_button");
    gtk_widget_set_tooltip_text(pt_button, "Periodic Table");
    gtk_grid_attach(GTK_GRID(education_grid), pt_button, 0, 2, 1, 1);
    g_signal_connect(pt_button, "clicked", G_CALLBACK(on_pt_button_clicked), notebook); 
   
    // SEMANTIC button
    GtkWidget *semantic_button = gtk_button_new_with_label("");
    gtk_widget_set_name(semantic_button, "semantic_button");
    gtk_widget_set_tooltip_text(semantic_button, "Semantic Scholar");
    gtk_grid_attach(GTK_GRID(education_grid), semantic_button, 1, 2, 1, 1);
    g_signal_connect(semantic_button, "clicked", G_CALLBACK(on_semantic_button_clicked), notebook);     
   
    // SCIENCEDIRECT button
    GtkWidget *sciencedirect_button = gtk_button_new_with_label("");
    gtk_widget_set_name(sciencedirect_button, "sciencedirect_button");
    gtk_widget_set_tooltip_text(sciencedirect_button, "Science Direct");
    gtk_grid_attach(GTK_GRID(education_grid), sciencedirect_button, 2, 2, 1, 1);
    g_signal_connect(sciencedirect_button, "clicked", G_CALLBACK(on_sciencedirect_button_clicked), notebook);
   
    // SCIENCE.GOV button
    GtkWidget *sciencegov_button = gtk_button_new_with_label("");
    gtk_widget_set_name(sciencegov_button, "sciencegov_button");
    gtk_widget_set_tooltip_text(sciencegov_button, "Science.gov");
    gtk_grid_attach(GTK_GRID(education_grid), sciencegov_button, 0, 3, 1, 1);
    g_signal_connect(sciencegov_button, "clicked", G_CALLBACK(on_sciencegov_button_clicked), notebook);   

    // WORLDWIDEHISTORY button
    GtkWidget *wwh_button = gtk_button_new_with_label("");
    gtk_widget_set_name(wwh_button, "wwh_button");
    gtk_widget_set_tooltip_text(wwh_button, "World Wide History");
    gtk_grid_attach(GTK_GRID(education_grid), wwh_button, 1, 3, 1, 1);
    g_signal_connect(wwh_button, "clicked", G_CALLBACK(on_wwh_button_clicked), notebook); 
   
    // WORLDWIDETELESCOPE button
    GtkWidget *wwt_button = gtk_button_new_with_label("");
    gtk_widget_set_name(wwt_button, "wwt_button");
    gtk_widget_set_tooltip_text(wwt_button, "World Wide Telescope");
    gtk_grid_attach(GTK_GRID(education_grid), wwt_button, 2, 3, 1, 1);
    g_signal_connect(wwt_button, "clicked", G_CALLBACK(on_wwt_button_clicked), notebook);     
//END GLOBE LINK

// NOTEBAR
GtkWidget *notebar = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_widget_set_size_request(notebar, 700, 1); 
gtk_widget_set_name(notebar, "notebar");
gtk_box_pack_start(GTK_BOX(hbox), notebar, FALSE, TRUE, 0);

// Create a toolbar
GtkWidget *notetoolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
gtk_widget_set_name(notetoolbar, "notetoolbar");
gtk_box_pack_start(GTK_BOX(notebar), notetoolbar, FALSE, FALSE, 0);

// Create a grid
GtkWidget *notegrid = gtk_grid_new();
gtk_widget_set_name(notegrid, "notegrid");
gtk_container_set_border_width(GTK_CONTAINER(notegrid), 10);
gtk_box_pack_start(GTK_BOX(notetoolbar), notegrid, TRUE, TRUE, 0);

// ADDTEXT button
GtkWidget *buttonaddtext = gtk_button_new_with_label("+");
gtk_widget_set_name(buttonaddtext, "buttonaddtext");
gtk_widget_set_tooltip_text(GTK_WIDGET(buttonaddtext), "Add tab"); 
gtk_grid_attach(GTK_GRID(notegrid), buttonaddtext, 1, 0, 1, 1);

// DELETETEXT button
GtkWidget *buttondeletetext = gtk_button_new_with_label("-");
gtk_widget_set_name(buttondeletetext, "buttondeletetext");
gtk_widget_set_tooltip_text(GTK_WIDGET(buttondeletetext), "Close tab"); 
gtk_grid_attach(GTK_GRID(notegrid), buttondeletetext, 2, 0, 1, 1);
  
// OPENTEXT button
GtkWidget *buttonopentext = gtk_button_new_with_label("Open");
gtk_widget_set_name(buttonopentext, "buttonopentext");
gtk_widget_set_tooltip_text(GTK_WIDGET(buttonopentext), "Open file"); 
gtk_grid_attach(GTK_GRID(notegrid), buttonopentext, 3, 0, 1, 1);

// SAVEASTEXT button
GtkWidget *buttonsavetext = gtk_button_new_with_label("Save as");
gtk_widget_set_name(buttonsavetext, "buttonsavetext");
gtk_widget_set_tooltip_text(GTK_WIDGET(buttonsavetext), "Save as"); 
gtk_grid_attach(GTK_GRID(notegrid), buttonsavetext, 4, 0, 1, 1);

// LEFTTEXT button
GtkWidget *buttonlefttext = gtk_button_new_with_label("");
gtk_widget_set_name(buttonlefttext, "buttonlefttext");
gtk_widget_set_tooltip_text(GTK_WIDGET(buttonlefttext), "Align left"); 
gtk_grid_attach(GTK_GRID(notegrid), buttonlefttext, 5, 0, 1, 1);

// CENTERTEXT button
GtkWidget *buttoncentertext = gtk_button_new_with_label("");
gtk_widget_set_name(buttoncentertext, "buttoncentertext");
gtk_widget_set_tooltip_text(GTK_WIDGET(buttoncentertext), "Align center"); 
gtk_grid_attach(GTK_GRID(notegrid), buttoncentertext, 6, 0, 1, 1);

// LIST button
GtkWidget *buttonlisttext = gtk_button_new_with_label("");
gtk_widget_set_name(buttonlisttext, "buttonlisttext");
gtk_widget_set_tooltip_text(GTK_WIDGET(buttonlisttext), "Align center"); 
gtk_grid_attach(GTK_GRID(notegrid), buttonlisttext, 7, 0, 1, 1);

// COMMENT button
GtkWidget *buttoncommenttext = gtk_button_new_with_label("");
gtk_widget_set_name(buttoncommenttext, "buttoncommenttext");
gtk_widget_set_tooltip_text(GTK_WIDGET(buttoncommenttext), "Comment"); 
gtk_grid_attach(GTK_GRID(notegrid), buttoncommenttext, 8, 0, 1, 1);

// NOTE NOTEBOOK
GtkNotebook *note_notebook = GTK_NOTEBOOK(gtk_notebook_new());
gtk_notebook_set_scrollable(note_notebook, TRUE); // Make the tabs scrollable
gtk_widget_set_name(GTK_WIDGET(note_notebook), "note_notebook"); // Set the notebook's name
gtk_box_pack_start(GTK_BOX(notebar), GTK_WIDGET(note_notebook), TRUE, TRUE, 0);

// Create a scrolled window
GtkWidget *scrolled_window = gtk_scrolled_window_new(NULL, NULL);
gtk_widget_set_name(scrolled_window, "scrolled_window");
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_window), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

// Create a text view widget
GtkWidget *text_view = gtk_text_view_new();
gtk_widget_set_name(text_view, "text_view");
gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(text_view), TRUE);
gtk_container_add(GTK_CONTAINER(scrolled_window), text_view); // Add text view to scrolled window
g_signal_connect(G_OBJECT(text_view), "paste-clipboard", G_CALLBACK(on_paste_clipboard), NULL);


// Add the scrolled window to the notebook
GtkWidget *tab_label1 = gtk_label_new("Tab");
gtk_notebook_append_page(note_notebook, scrolled_window, tab_label1);

// Set maximum width to accommodate approximately 48 characters
int char_width = 12; // Adjust this based on your font
int max_width = char_width * 48; // 100 characters
gtk_widget_set_size_request(text_view, max_width, -1); // -1 means no limit on height

// Set maximum height to 400 pixels
gtk_widget_set_size_request(scrolled_window, 1, 900); // -1 means no limit on width

// Enable text wrapping
gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view), GTK_WRAP_WORD);

// Set some sample text
GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
gtk_text_buffer_set_text(buffer, "", -1);

// Show all widgets
gtk_widget_show_all(notebar);

//NOTEBBOK FUNCTION
g_signal_connect(buttonaddtext, "clicked", G_CALLBACK(create_new_tab), note_notebook);    
g_signal_connect(buttondeletetext, "clicked", G_CALLBACK(close_active_notetab), note_notebook); 
g_signal_connect(buttonopentext, "clicked", G_CALLBACK(on_buttonopentext_clicked), note_notebook);  
g_signal_connect(buttonsavetext, "clicked", G_CALLBACK(on_buttonsavetext_clicked), note_notebook);
g_signal_connect(buttonlefttext, "clicked", G_CALLBACK(on_left_align_text), note_notebook);
g_signal_connect(buttoncentertext, "clicked", G_CALLBACK(on_center_text), note_notebook);
g_signal_connect(buttonlisttext, "clicked", G_CALLBACK(on_list_text), note_notebook);
g_signal_connect(buttoncommenttext, "clicked", G_CALLBACK(on_comment_text), note_notebook);
g_signal_connect(text_view, "key-press-event", G_CALLBACK(on_enter_pressed), note_notebook);

//NOTEBOOK END

//////////BOOKMARKS
GtkWidget *bookmarkbar = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_widget_set_size_request(bookmarkbar, 240, 1); 
gtk_widget_set_name(bookmarkbar, "bookmarkbar");
gtk_box_pack_start(GTK_BOX(hbox), bookmarkbar, FALSE, TRUE, 0);

GtkWidget *bookmarkbar_title = gtk_label_new("BOOKMARKS");
gtk_widget_set_name(bookmarkbar_title, "downloadbar_title");
gtk_box_pack_start(GTK_BOX(bookmarkbar), bookmarkbar_title, FALSE, FALSE, 0);

GtkWidget *addbookmark_button = gtk_button_new_with_label("ADD");
gtk_widget_set_name(addbookmark_button, "add_button");
gtk_widget_set_halign(addbookmark_button, GTK_ALIGN_CENTER);
gtk_widget_set_valign(addbookmark_button, GTK_ALIGN_CENTER);

gtk_box_pack_start(GTK_BOX(bookmarkbar), addbookmark_button, FALSE, FALSE, 0);
gtk_widget_show(addbookmark_button);

    // Set user data to be accessed in callback
    g_object_set_data(G_OBJECT(addbookmark_button), "bookmarkbar", bookmarkbar);

    // Connect the add button signal to the callback function
    g_signal_connect(addbookmark_button, "clicked", G_CALLBACK(add_bookmark), notebook);

// Set the notebook data on the bookmarkbar
set_notebook_data(bookmarkbar, notebook);



//DOWNLOAD BAR
// Create a vertical box for the menu
GtkWidget *downloadbar = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_widget_set_size_request(downloadbar, 450, 1); // Set width to 48 pixels
gtk_widget_set_name(downloadbar, "downloadbar");
gtk_box_pack_start(GTK_BOX(hbox), downloadbar, FALSE, TRUE, 0);


GtkWidget *downloadbar_title = gtk_label_new("DOWNLOAD");
gtk_widget_set_name(downloadbar_title, "downloadbar_title");
gtk_box_pack_start(GTK_BOX(downloadbar), downloadbar_title, FALSE, FALSE, 0);
/*
WebKitWebContext *context = webkit_web_context_get_default();
g_signal_connect(context, "download-started", G_CALLBACK(download_started_callback), NULL);*/
//g_print("Connected download-started signal on WebKitWebContext: %p\n", (void *)context);
g_signal_connect(webkit_web_view_get_context(web_view), "download-started", G_CALLBACK(download_started_callback), NULL);

GtkWidget *opendownloadfolder_button = gtk_button_new_with_label("Open Download Folder");
gtk_widget_set_name(opendownloadfolder_button, "add_button");
gtk_widget_set_size_request(opendownloadfolder_button, 150, -1);
gtk_widget_set_hexpand(opendownloadfolder_button, FALSE);
gtk_widget_set_halign(opendownloadfolder_button, GTK_ALIGN_CENTER);
gtk_box_pack_start(GTK_BOX(downloadbar), opendownloadfolder_button, FALSE, FALSE, 0);
g_signal_connect(opendownloadfolder_button, "button-press-event", G_CALLBACK(open_download_folder), NULL);

download_event_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_box_pack_start(GTK_BOX(downloadbar), download_event_box, FALSE, FALSE, 0);

// Create a scrolled window
GtkWidget *download_scrolled_window = gtk_scrolled_window_new(NULL, NULL);
gtk_widget_set_name(download_scrolled_window, "scrolled_window");
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(download_scrolled_window), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(downloadbar), download_scrolled_window, TRUE, TRUE, 0); 

// Create a list box to hold the file list
file_list = gtk_list_box_new();
gtk_widget_set_name(file_list, "file_list");
// Add the file list to the scrolled window
gtk_container_add(GTK_CONTAINER(download_scrolled_window), file_list);
// Initial display of the Downloads folder contents
show_downloads(file_list);
// Retrieve the top-level widget after adding to the container
GtkWidget *toplevel = gtk_widget_get_toplevel(file_list);



// HISTORY BAR
GtkWidget *historybar = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_widget_set_size_request(historybar, 400, 1);
gtk_widget_set_name(historybar, "historybar");
gtk_box_pack_start(GTK_BOX(hbox), historybar, false, true, 0);

GtkWidget *web_history_label = gtk_label_new("WEB HISTORY");
gtk_widget_set_halign(web_history_label, GTK_ALIGN_CENTER);
gtk_widget_set_name(web_history_label, "web_history_label");
gtk_box_pack_start(GTK_BOX(historybar), web_history_label, false, false, 0);

GtkWidget *history_scrolled_window = gtk_scrolled_window_new(NULL, NULL);
gtk_widget_set_name(history_scrolled_window, "scrolled_window");
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(history_scrolled_window), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

GtkWidget *history_container = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_container_add(GTK_CONTAINER(history_scrolled_window), history_container);

gtk_box_pack_start(GTK_BOX(historybar), history_scrolled_window, TRUE, TRUE, 0);

set_history_notebook_data(history_container, notebook);
load_history(history_container);

gtk_widget_show_all(historybar);

//CALCULATOR
GtkWidget *calculatorbar = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_widget_set_size_request(calculatorbar, 500, 1);
gtk_widget_set_hexpand(calculatorbar, FALSE);
gtk_widget_set_name(calculatorbar, "calculatorbar");
gtk_box_pack_start(GTK_BOX(hbox), calculatorbar, false, true, 0);

GtkWidget *calculator_app_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
gtk_widget_set_size_request(calculator_app_box, 500, 740);
gtk_widget_set_hexpand(calculator_app_box, FALSE);
gtk_widget_set_vexpand(calculator_app_box, FALSE);
gtk_widget_set_name(calculator_app_box, "calculatorbar");
gtk_box_pack_start(GTK_BOX(calculatorbar), calculator_app_box, false, false, 0);



    // Create a vertical box
    GtkWidget *calculatorbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
    gtk_container_add(GTK_CONTAINER(calculator_app_box), calculatorbox);


    // Call the setup function
    setup_drawing_area(calculatorbox);

    // Create a grid
    GtkWidget *calcgrid = gtk_grid_new();
    gtk_widget_set_hexpand(calcgrid, TRUE); 
    gtk_grid_set_column_homogeneous(GTK_GRID(calcgrid), TRUE);
    gtk_box_pack_start(GTK_BOX(calculatorbox), calcgrid, TRUE, TRUE, 0);

    // Create entry widgets
    GtkWidget *entry_expression = gtk_entry_new();
    gtk_widget_set_name(entry_expression, "entry_expression");
    gtk_widget_set_hexpand(entry_expression, TRUE); 
    gtk_grid_attach(GTK_GRID(calcgrid), entry_expression, 0, 0, 7, 1);

    GtkWidget *entry_result = gtk_entry_new();
    gtk_widget_set_name(entry_result, "entry_result");
    gtk_widget_set_hexpand(entry_result, TRUE); 
    gtk_grid_attach(GTK_GRID(calcgrid), entry_result, 0, 1, 7, 1);


// CALCULATOR BUTTONS
GtkWidget *buttongrid = gtk_grid_new();
gtk_container_set_border_width(GTK_CONTAINER(buttongrid), 13); // Set border width for buttongrid
gtk_grid_set_column_homogeneous(GTK_GRID(buttongrid), TRUE);

// FIRST ROW
GtkWidget *squareroot = gtk_button_new_with_label("²√");
gtk_widget_set_name(squareroot, "calcbuttons");
gtk_grid_attach(GTK_GRID(buttongrid), squareroot, 0, 0, 1, 1); // Column 0, Row 0
g_signal_connect(squareroot, "clicked", G_CALLBACK(on_squareroot_button_clicked), entry_expression);
g_object_set_data(G_OBJECT(squareroot), "entry_result", entry_result); // Pass entry_result to the square root function

GtkWidget *cuberoot = gtk_button_new_with_label("³√");
gtk_widget_set_name(cuberoot, "calcbuttons");
gtk_grid_attach(GTK_GRID(buttongrid), cuberoot, 1, 0, 1, 1); // Column 1, Row 0
g_signal_connect(cuberoot, "clicked", G_CALLBACK(on_cuberoot_button_clicked), entry_expression);
g_object_set_data(G_OBJECT(cuberoot), "entry_result", entry_result); // Pass entry_result to the square root function

GtkWidget *xroot = gtk_button_new_with_label("^√");
gtk_widget_set_name(xroot, "calcbuttons");
gtk_grid_attach(GTK_GRID(buttongrid), xroot, 2, 0, 1, 1); // Column 0, Row 0
g_signal_connect(xroot, "clicked", G_CALLBACK(on_xroot_button_clicked), entry_expression);
g_object_set_data(G_OBJECT(xroot), "entry_result", entry_result); // Pass entry_result to the xroot function



GtkWidget *ac = gtk_button_new_with_label("AC");
gtk_widget_set_name(ac, "ac");
gtk_grid_attach(GTK_GRID(buttongrid), ac, 3, 0, 1, 1); // Column 1, Row 0
g_signal_connect(ac, "clicked", G_CALLBACK(on_ac_button_clicked), entry_expression);
g_object_set_data(G_OBJECT(ac), "entry_result", entry_result); // Pass entry_result to the AC function
g_object_set_data(G_OBJECT(ac), "calcdrawing_area", calcdrawing_area); // Pass drawing_area to the AC function


GtkWidget *pgreco = gtk_button_new_with_label("π");
gtk_widget_set_name(pgreco, "calcbuttons");
g_signal_connect(pgreco, "clicked", G_CALLBACK(on_pgreco_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), pgreco, 4, 0, 1, 1); // Column 1, Row 0

GtkWidget *eulero = gtk_button_new_with_label("e");
gtk_widget_set_name(eulero, "calcbuttons");
g_signal_connect(eulero, "clicked", G_CALLBACK(on_eulero_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), eulero, 5, 0, 1, 1); // Column 0, Row 056

GtkWidget *cosh = gtk_button_new_with_label("cosh");
gtk_widget_set_name(cosh, "calcbuttons");
g_signal_connect(cosh, "clicked", G_CALLBACK(on_cosh_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), cosh, 6, 0, 1, 1); // Column 1, Row 0

//SECOND ROW //    "7", "8", "9", "/", "x^", "log", "sinh",
GtkWidget *seven = gtk_button_new_with_label("7");
gtk_widget_set_name(seven, "calcbuttons");
g_signal_connect(seven, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), seven, 0, 1, 1, 1); // Column 0, Row 0

GtkWidget *eight = gtk_button_new_with_label("8");
gtk_widget_set_name(eight, "calcbuttons");
g_signal_connect(eight, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), eight, 1, 1, 1, 1); // Column 1, Row 0

GtkWidget *nine = gtk_button_new_with_label("9");
gtk_widget_set_name(nine, "calcbuttons");
g_signal_connect(nine, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), nine, 2, 1, 1, 1); // Column 0, Row 0

GtkWidget *division = gtk_button_new_with_label("/");
gtk_widget_set_name(division, "calcbuttons");
gtk_grid_attach(GTK_GRID(buttongrid), division, 3, 1, 1, 1); // Column 1, Row 0
g_signal_connect(division, "clicked", G_CALLBACK(on_division_button_clicked), entry_expression);


GtkWidget *xpower = gtk_button_new_with_label("x^");
gtk_widget_set_name(xpower, "calcbuttons");
g_signal_connect(xpower, "clicked", G_CALLBACK(on_xpower_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), xpower, 4, 1, 1, 1); // Column 1, Row 0

GtkWidget *log = gtk_button_new_with_label("log");
gtk_widget_set_name(log, "calcbuttons");
g_signal_connect(log, "clicked", G_CALLBACK(on_log_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), log, 5, 1, 1, 1); // Column 0, Row 0

GtkWidget *sinh = gtk_button_new_with_label("sinh");
gtk_widget_set_name(sinh, "calcbuttons");
g_signal_connect(sinh , "clicked", G_CALLBACK(on_sinh_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), sinh, 6, 1, 1, 1); // Column 1, Row 0

//THIRD ROW //     "4", "5", "6", "*", "x³", "cos", "acos",
GtkWidget *four = gtk_button_new_with_label("4");
gtk_widget_set_name(four, "calcbuttons");
g_signal_connect(four, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), four, 0, 2, 1, 1); // Column 0, Row 0

GtkWidget *five = gtk_button_new_with_label("5");
gtk_widget_set_name(five, "calcbuttons");
g_signal_connect(five, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), five, 1, 2, 1, 1); // Column 1, Row 0

GtkWidget *six = gtk_button_new_with_label("6");
gtk_widget_set_name(six, "calcbuttons");
g_signal_connect(six, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), six, 2, 2, 1, 1); // Column 0, Row 0

GtkWidget *moltiplication = gtk_button_new_with_label("*");
gtk_widget_set_name(moltiplication, "calcbuttons");
gtk_grid_attach(GTK_GRID(buttongrid), moltiplication, 3, 2, 1, 1); // Column 1, Row 0
g_signal_connect(moltiplication, "clicked", G_CALLBACK(on_multiplication_button_clicked), entry_expression);


GtkWidget *cubepower = gtk_button_new_with_label("x³");
gtk_widget_set_name(cubepower, "calcbuttons");
g_signal_connect(cubepower , "clicked", G_CALLBACK(on_cubepower_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), cubepower, 4, 2, 1, 1); // Column 1, Row 0

GtkWidget *cos = gtk_button_new_with_label("cos");
gtk_widget_set_name(cos, "calcbuttons");
g_signal_connect(cos, "clicked", G_CALLBACK(on_cos_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), cos, 5, 2, 1, 1); // Column 0, Row 0

GtkWidget *acos = gtk_button_new_with_label("acos");
gtk_widget_set_name(acos, "calcbuttons");
g_signal_connect(acos, "clicked", G_CALLBACK(on_acos_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), acos, 6, 2, 1, 1); // Column 1, Row 0

//FOURTH ROW //      "1", "2", "3", "-", "x²", "sin", "asin",
GtkWidget *one = gtk_button_new_with_label("1");
gtk_widget_set_name(one, "calcbuttons");
g_signal_connect(one, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), one, 0, 3, 1, 1); // Column 0, Row 0

GtkWidget *two = gtk_button_new_with_label("2");
gtk_widget_set_name(two, "calcbuttons");
g_signal_connect(two, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), two, 1, 3, 1, 1); // Column 1, Row 0

GtkWidget *three = gtk_button_new_with_label("3");
gtk_widget_set_name(three, "calcbuttons");
g_signal_connect(three, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), three, 2, 3, 1, 1); // Column 0, Row 0

GtkWidget *subtraction = gtk_button_new_with_label("-");
gtk_widget_set_name(subtraction, "calcbuttons");
g_signal_connect(subtraction, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
g_signal_connect(subtraction, "clicked", G_CALLBACK(on_subtraction_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), subtraction, 3, 3, 1, 1); // Column 1, Row 0

GtkWidget *squarepower = gtk_button_new_with_label("x²");
gtk_widget_set_name(squarepower, "calcbuttons");
g_signal_connect(squarepower , "clicked", G_CALLBACK(on_squarepower_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), squarepower, 4, 3, 1, 1); // Column 1, Row 0

GtkWidget *sin = gtk_button_new_with_label("sin");
gtk_widget_set_name(sin, "calcbuttons");
g_signal_connect(sin , "clicked", G_CALLBACK(on_sin_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), sin, 5, 3, 1, 1); // Column 0, Row 0

GtkWidget *asin = gtk_button_new_with_label("asin");
gtk_widget_set_name(asin, "calcbuttons");
g_signal_connect(asin , "clicked", G_CALLBACK(on_asin_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), asin, 6, 3, 1, 1); // Column 1, Row 0

//FIFTH ROW //       "0", ".", "=", "+", "x-¹", "tan", "atan",
GtkWidget *zero = gtk_button_new_with_label("0");
gtk_widget_set_name(zero, "calcbuttons");
g_signal_connect(zero, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), zero, 0, 4, 1, 1); // Column 0, Row 0

GtkWidget *point = gtk_button_new_with_label(".");
gtk_widget_set_name(point, "calcbuttons");
g_signal_connect(point, "clicked", G_CALLBACK(on_number_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), point, 1, 4, 1, 1); // Column 1, Row 0

GtkWidget *equal = gtk_button_new_with_label("=");
gtk_widget_set_name(equal, "calcbuttons");
gtk_grid_attach(GTK_GRID(buttongrid), equal, 2, 4, 1, 1); // Column 0, Row 0

GtkWidget *plus = gtk_button_new_with_label("+");
gtk_widget_set_name(plus, "calcbuttons");
gtk_grid_attach(GTK_GRID(buttongrid), plus, 3, 4, 1, 1); // Column 1, Row 0
g_signal_connect(plus, "clicked", G_CALLBACK(on_plus_button_clicked), entry_expression);


GtkWidget *negativepower = gtk_button_new_with_label("x-¹");
gtk_widget_set_name(negativepower, "calcbuttons");
g_signal_connect(negativepower, "clicked", G_CALLBACK(on_negativepower_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), negativepower, 4, 4, 1, 1); // Column 1, Row 0

GtkWidget *tan = gtk_button_new_with_label("tan");
gtk_widget_set_name(tan, "calcbuttons");
g_signal_connect(tan, "clicked", G_CALLBACK(on_tan_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), tan, 5, 4, 1, 1); // Column 0, Row 0

GtkWidget *atan = gtk_button_new_with_label("atan");
gtk_widget_set_name(atan, "calcbuttons");
g_signal_connect(atan, "clicked", G_CALLBACK(on_atan_button_clicked), entry_expression);
gtk_grid_attach(GTK_GRID(buttongrid), atan, 6, 4, 1, 1); // Column 1, Row 0

// Pack buttongrid into calculatorbox
gtk_box_pack_start(GTK_BOX(calculatorbox), buttongrid, TRUE, TRUE, 0);

//MORE START
GtkWidget *morescrolled_window = gtk_scrolled_window_new(NULL, NULL);
gtk_box_pack_start(GTK_BOX(calculator_app_box), morescrolled_window, FALSE, FALSE, 0);
gtk_widget_set_size_request(morescrolled_window, 220, 740);
    // MOREBOX
    morebox =gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
    gtk_container_add(GTK_CONTAINER(morescrolled_window), morebox);

/*GtkWidget *calcmorebutton = gtk_button_new_with_label("MORE");
gtk_widget_set_name(calcmorebutton, "calcmorebutton");
gtk_grid_attach(GTK_GRID(buttongrid), calcmorebutton, 1, 7, 5, 1); // Column 1, Row 0
g_signal_connect(calcmorebutton, "clicked", G_CALLBACK(on_calcmorebutton_button_clicked), morescrolled_window);
*/



more_title = gtk_label_new("MORE");
gtk_widget_set_name(more_title, "more_title");
gtk_box_pack_start(GTK_BOX(morebox), more_title, FALSE, FALSE, 0);

//LINEAR EQUATION
GtkWidget *line_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
gtk_widget_set_name(line_box, "moreboxes");
gtk_container_add(GTK_CONTAINER(morebox), line_box);

// Create the label and entries
label_line = gtk_label_new("Linear equation");
gtk_widget_set_name(label_line, "more_label_title");
gtk_box_pack_start(GTK_BOX(line_box), label_line, FALSE, FALSE, 0);

label_line_formula = gtk_label_new("y = mx + b");
gtk_widget_set_name(label_line_formula, "more_label_title");
gtk_box_pack_start(GTK_BOX(line_box), label_line_formula, FALSE, FALSE, 0);

label_m = gtk_label_new("Enter the value for m:");
gtk_widget_set_name(label_m, "enter_label");
gtk_box_pack_start(GTK_BOX(line_box), label_m, FALSE, FALSE, 0);

entry_m = gtk_entry_new();
gtk_widget_set_name(entry_m, "more_entry");
gtk_box_pack_start(GTK_BOX(line_box), entry_m, FALSE, FALSE, 0);

label_b = gtk_label_new("Enter the value for b:");
gtk_widget_set_name(label_b, "enter_label");
gtk_box_pack_start(GTK_BOX(line_box), label_b, FALSE, FALSE, 0);

entry_b = gtk_entry_new();
gtk_widget_set_name(entry_b, "more_entry");
gtk_box_pack_start(GTK_BOX(line_box), entry_b, FALSE, FALSE, 0);

// Store the entries in an array
GtkWidget *entries_line[2] = {entry_m, entry_b};

// RESULT
label_result_line = gtk_label_new("");
gtk_widget_set_name(label_result_line, "result_label");
gtk_box_pack_start(GTK_BOX(line_box), label_result_line, FALSE, FALSE, 0);

// M
label_m_line = gtk_label_new("");
gtk_widget_set_name(label_m_line, "result_label");
gtk_box_pack_start(GTK_BOX(line_box), label_m_line, FALSE, FALSE, 0);

// B
label_b_line = gtk_label_new("");
gtk_widget_set_name(label_b_line, "result_label");
gtk_box_pack_start(GTK_BOX(line_box), label_b_line, FALSE, FALSE, 0);

gpointer data_line[4] = {entries_line, label_m_line, label_b_line, label_result_line};

GtkWidget *show_line = gtk_button_new_with_label("Show line");
gtk_widget_set_name(show_line, "calcmorebutton");
g_signal_connect(show_line, "clicked", G_CALLBACK(on_show_line_button_clicked), data_line);
gtk_box_pack_start(GTK_BOX(line_box), show_line, FALSE, FALSE, 0);

// Create and connect the "linear equation" button
GtkWidget *line = gtk_button_new_with_label("linear equation");
gtk_widget_set_name(line, "morebutton");
g_signal_connect(line, "clicked", G_CALLBACK(on_line_button_clicked), line_box);
gtk_box_pack_start(GTK_BOX(morebox), line, FALSE, FALSE, 0);
///LINEAR EQUATION END


    // PARABOLA Y
    GtkWidget *parabolabox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
    gtk_widget_set_name(parabolabox, "moreboxes");
    gtk_container_add(GTK_CONTAINER(morebox), parabolabox);

    // Create the label and entries
    labelparabola = gtk_label_new("parabola Y");
    gtk_widget_set_name(labelparabola, "more_label_title");
    gtk_box_pack_start(GTK_BOX(parabolabox), labelparabola, FALSE, FALSE, 0);

    labelparabola_y_formula = gtk_label_new("y = ax² + bx + c");
    gtk_widget_set_name(labelparabola_y_formula, "more_label_title");
    gtk_box_pack_start(GTK_BOX(parabolabox), labelparabola_y_formula, FALSE, FALSE, 0);

    label_a_parabola = gtk_label_new("Enter the value for a:");
    gtk_widget_set_name(label_a_parabola, "enter_label");
    gtk_box_pack_start(GTK_BOX(parabolabox), label_a_parabola, FALSE, FALSE, 0);

    entrya_parabolay = gtk_entry_new();
    gtk_widget_set_name(entrya_parabolay, "more_entry");
    gtk_box_pack_start(GTK_BOX(parabolabox), entrya_parabolay, FALSE, FALSE, 0);

    label_b_parabola = gtk_label_new("Enter the value for b:");
    gtk_widget_set_name(label_b_parabola, "enter_label");
    gtk_box_pack_start(GTK_BOX(parabolabox), label_b_parabola, FALSE, FALSE, 0);

    entryb_parabolay = gtk_entry_new();
    gtk_widget_set_name(entryb_parabolay, "more_entry");
    gtk_box_pack_start(GTK_BOX(parabolabox), entryb_parabolay, FALSE, FALSE, 0);

    label_c_parabola = gtk_label_new("Enter the value for c:");
    gtk_widget_set_name(label_c_parabola, "enter_label");
    gtk_box_pack_start(GTK_BOX(parabolabox), label_c_parabola, FALSE, FALSE, 0);

    entryc_parabolay = gtk_entry_new();
    gtk_widget_set_name(entryc_parabolay, "more_entry");
    gtk_box_pack_start(GTK_BOX(parabolabox), entryc_parabolay, FALSE, FALSE, 0);

   // Store the entries in an array
    GtkWidget *entries[3] = {entrya_parabolay, entryb_parabolay, entryc_parabolay};
    
//RISULTATO
    label_result_parabolay = gtk_label_new("");
    gtk_widget_set_name( label_result_parabolay, "result_label");
    gtk_box_pack_start(GTK_BOX(parabolabox),  label_result_parabolay, FALSE, FALSE, 0);
//DELTA
    label_delta_parabolay = gtk_label_new("");
    gtk_widget_set_name(  label_delta_parabolay, "result_label");
    gtk_box_pack_start(GTK_BOX(parabolabox),  label_delta_parabolay, FALSE, FALSE, 0);
//X VERTEX
    label_vertex_x1_parabolay = gtk_label_new("");
    gtk_widget_set_name(  label_vertex_x1_parabolay, "result_label");
    gtk_box_pack_start(GTK_BOX(parabolabox),  label_vertex_x1_parabolay, FALSE, FALSE, 0);
    
    label_vertex_x2_parabolay = gtk_label_new("");
    gtk_widget_set_name(  label_vertex_x2_parabolay, "result_label");
    gtk_box_pack_start(GTK_BOX(parabolabox),  label_vertex_x2_parabolay, FALSE, FALSE, 0);
//Y VERTEX
    label_vertex_y_parabolay = gtk_label_new("");
    gtk_widget_set_name(  label_vertex_y_parabolay, "result_label");
    gtk_box_pack_start(GTK_BOX(parabolabox),  label_vertex_y_parabolay, FALSE, FALSE, 0);
//X FOCUS
    label_focus_x_parabolay = gtk_label_new("");
    gtk_widget_set_name(  label_focus_x_parabolay, "result_label");
    gtk_box_pack_start(GTK_BOX(parabolabox),  label_focus_x_parabolay, FALSE, FALSE, 0);
//Y FOCUS
    label_focus_y_parabolay = gtk_label_new("");
    gtk_widget_set_name(  label_focus_y_parabolay, "result_label");
    gtk_box_pack_start(GTK_BOX(parabolabox),  label_focus_y_parabolay, FALSE, FALSE, 0);
//AXYS OF SYMMETRY
    sim_axys_parabolay = gtk_label_new("");
    gtk_widget_set_name(  sim_axys_parabolay, "result_label");
    gtk_box_pack_start(GTK_BOX(parabolabox),  sim_axys_parabolay, FALSE, FALSE, 0);
//DIRECTRIX
    directrix = gtk_label_new("");
    gtk_widget_set_name(  directrix, "result_label");
    gtk_box_pack_start(GTK_BOX(parabolabox),  directrix, FALSE, FALSE, 0);
    
    gpointer data[10] = {entries, label_result_parabolay, label_vertex_x1_parabolay, label_vertex_x2_parabolay, label_delta_parabolay, label_vertex_y_parabolay, label_focus_x_parabolay, label_focus_y_parabolay, sim_axys_parabolay, directrix};
   
   // Create and connect the "Show parabola" button
    GtkWidget *show_parabolay = gtk_button_new_with_label("Show parabola");
    gtk_widget_set_name(show_parabolay, "calcmorebutton");
    g_signal_connect(show_parabolay, "clicked", G_CALLBACK(on_show_parabolay_button_clicked), data);
    gtk_box_pack_start(GTK_BOX(parabolabox), show_parabolay, FALSE, FALSE, 0);


    // Create and connect the "parabola" button
    GtkWidget *parabola = gtk_button_new_with_label("parabola y");
    gtk_widget_set_name(parabola, "morebutton");
    g_signal_connect(parabola, "clicked", G_CALLBACK(on_parabola_button_clicked), parabolabox);
    gtk_box_pack_start(GTK_BOX(morebox), parabola, FALSE, FALSE, 0);
////END PARABOLA Y

///////////PARABOLA X
    GtkWidget *parabola_x_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
    gtk_widget_set_name(parabola_x_box, "moreboxes");
    gtk_container_add(GTK_CONTAINER(morebox), parabola_x_box);

    // Create the label and entries
    labelparabolax = gtk_label_new("parabola X");
    gtk_widget_set_name(labelparabolax, "more_label_title");
    gtk_box_pack_start(GTK_BOX(parabola_x_box), labelparabolax, FALSE, FALSE, 0);

    labelparabola_x_formula = gtk_label_new("x = ay² + by + c");
    gtk_widget_set_name(labelparabola_x_formula, "more_label_title");
    gtk_box_pack_start(GTK_BOX(parabola_x_box), labelparabola_x_formula, FALSE, FALSE, 0);
 
    label_a_parabola = gtk_label_new("Enter the value for a:");
    gtk_widget_set_name(label_a_parabola, "enter_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box), label_a_parabola, FALSE, FALSE, 0);
  
    entrya_parabolax = gtk_entry_new();
    gtk_widget_set_name(entrya_parabolax, "more_entry");
    gtk_box_pack_start(GTK_BOX(parabola_x_box), entrya_parabolax, FALSE, FALSE, 0);
    
    label_b_parabola = gtk_label_new("Enter the value for b:");
    gtk_widget_set_name(label_b_parabola, "enter_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box), label_b_parabola, FALSE, FALSE, 0);

    entryb_parabolax = gtk_entry_new();
    gtk_widget_set_name(entryb_parabolax, "more_entry");
    gtk_box_pack_start(GTK_BOX(parabola_x_box), entryb_parabolax, FALSE, FALSE, 0);

    label_c_parabola = gtk_label_new("Enter the value for c:");
    gtk_widget_set_name(label_c_parabola, "enter_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box), label_c_parabola, FALSE, FALSE, 0);

    entryc_parabolax = gtk_entry_new();
    gtk_widget_set_name(entryc_parabolax, "more_entry");
    gtk_box_pack_start(GTK_BOX(parabola_x_box), entryc_parabolax, FALSE, FALSE, 0);
 
   // Store the entries in an array
    GtkWidget *entriesx[3] = {entrya_parabolax, entryb_parabolax, entryc_parabolax};

//RISULTATO

    label_result_parabolax = gtk_label_new("");
    gtk_widget_set_name(label_result_parabolax, "result_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box), label_result_parabolax, FALSE, FALSE, 0);
//DELTA
    label_delta_parabolax = gtk_label_new("");
    gtk_widget_set_name(  label_delta_parabolax, "result_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box),  label_delta_parabolax, FALSE, FALSE, 0);
//X VERTEX
    label_vertex_x_parabolax = gtk_label_new("");
    gtk_widget_set_name(  label_vertex_x_parabolax, "result_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box),  label_vertex_x_parabolax, FALSE, FALSE, 0);
//Y VERTEX
    label_vertex_y1_parabolax = gtk_label_new("");
    gtk_widget_set_name(  label_vertex_y1_parabolax, "result_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box),  label_vertex_y1_parabolax, FALSE, FALSE, 0);
    
    label_vertex_y2_parabolax = gtk_label_new("");
    gtk_widget_set_name(  label_vertex_y2_parabolax, "result_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box),  label_vertex_y2_parabolax, FALSE, FALSE, 0);    
//X FOCUS
    label_focus_x_parabolax = gtk_label_new("");
    gtk_widget_set_name(  label_focus_x_parabolax, "result_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box),  label_focus_x_parabolax, FALSE, FALSE, 0);
//Y FOCUS
    label_focus_y_parabolax = gtk_label_new("");
    gtk_widget_set_name(  label_focus_y_parabolax, "result_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box),  label_focus_y_parabolax, FALSE, FALSE, 0);
//AXYS OF SYMMETRY
    sim_axys_parabolax = gtk_label_new("");
    gtk_widget_set_name(  sim_axys_parabolax, "result_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box),  sim_axys_parabolax, FALSE, FALSE, 0);
//DIRECTRIX
    directrix_x = gtk_label_new("");
    gtk_widget_set_name(  directrix_x, "result_label");
    gtk_box_pack_start(GTK_BOX(parabola_x_box),  directrix_x, FALSE, FALSE, 0);
    
gpointer datax[10] = {entriesx, label_result_parabolax, label_vertex_x_parabolax, label_delta_parabolax, label_vertex_y1_parabolax, label_vertex_y2_parabolax, label_focus_x_parabolax, label_focus_y_parabolax, sim_axys_parabolax, directrix_x};

 // Create and connect the "Show parabola" button
    GtkWidget *show_parabolax = gtk_button_new_with_label("Show parabola");
    gtk_widget_set_name(show_parabolax, "calcmorebutton");
    g_signal_connect(show_parabolax, "clicked", G_CALLBACK(on_show_parabolax_button_clicked), datax);
    gtk_box_pack_start(GTK_BOX(parabola_x_box), show_parabolax, FALSE, FALSE, 0);


    // Create and connect the "parabola" button
    GtkWidget *parabolax = gtk_button_new_with_label("parabola x");
    gtk_widget_set_name(parabolax, "morebutton");
    g_signal_connect(parabolax, "clicked", G_CALLBACK(on_parabolax_button_clicked), parabola_x_box);
    gtk_box_pack_start(GTK_BOX(morebox), parabolax, FALSE, FALSE, 0);
///END PARABOLA X

//SQUARE EQUATION
GtkWidget *eq_second_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
gtk_widget_set_name(eq_second_box, "moreboxes");
gtk_container_add(GTK_CONTAINER(morebox), eq_second_box);

// Create the label and entries
GtkWidget *label_eq_second = gtk_label_new("Square function");
gtk_widget_set_name(label_eq_second, "more_label_title");
gtk_box_pack_start(GTK_BOX(eq_second_box), label_eq_second, FALSE, FALSE, 0);

GtkWidget *label_eq_second_formula = gtk_label_new("ax² + bx + c = 0");
gtk_widget_set_name(label_eq_second_formula, "more_label_title");
gtk_box_pack_start(GTK_BOX(eq_second_box), label_eq_second_formula, FALSE, FALSE, 0);

GtkWidget *label_a_eqsecond = gtk_label_new("Enter the value for a:");
gtk_widget_set_name(label_a_eqsecond, "enter_label");
gtk_box_pack_start(GTK_BOX(eq_second_box), label_a_eqsecond, FALSE, FALSE, 0);

GtkWidget *entrya_eqsecond = gtk_entry_new();
gtk_widget_set_name(entrya_eqsecond, "more_entry");
gtk_box_pack_start(GTK_BOX(eq_second_box), entrya_eqsecond, FALSE, FALSE, 0);

GtkWidget *label_b_eqsecond = gtk_label_new("Enter the value for b:");
gtk_widget_set_name(label_b_eqsecond, "enter_label");
gtk_box_pack_start(GTK_BOX(eq_second_box), label_b_eqsecond, FALSE, FALSE, 0);

GtkWidget *entryb_eqsecond = gtk_entry_new();
gtk_widget_set_name(entryb_eqsecond, "more_entry");
gtk_box_pack_start(GTK_BOX(eq_second_box), entryb_eqsecond, FALSE, FALSE, 0);

GtkWidget *label_c_eqsecond = gtk_label_new("Enter the value for c:");
gtk_widget_set_name(label_c_eqsecond, "enter_label");
gtk_box_pack_start(GTK_BOX(eq_second_box), label_c_eqsecond, FALSE, FALSE, 0);

GtkWidget *entryc_eqsecond = gtk_entry_new();
gtk_widget_set_name(entryc_eqsecond, "more_entry");
gtk_box_pack_start(GTK_BOX(eq_second_box), entryc_eqsecond, FALSE, FALSE, 0);
    
// Store the entries in an array
GtkWidget *entries_eq2[3] = {entrya_eqsecond, entryb_eqsecond, entryc_eqsecond};

GtkWidget *label_result_eqsecond = gtk_label_new("");
gtk_widget_set_name(label_result_eqsecond, "result_label");
gtk_box_pack_start(GTK_BOX(eq_second_box), label_result_eqsecond, FALSE, FALSE, 0);

GtkWidget *label_delta_eqsecond = gtk_label_new("");
gtk_widget_set_name(label_delta_eqsecond, "result_label");
gtk_box_pack_start(GTK_BOX(eq_second_box), label_delta_eqsecond, FALSE, FALSE, 0);

GtkWidget *label_x1_eqsecond = gtk_label_new("");
gtk_widget_set_name(label_x1_eqsecond, "result_label");
gtk_box_pack_start(GTK_BOX(eq_second_box), label_x1_eqsecond, FALSE, FALSE, 0);

GtkWidget *label_x2_eqsecond = gtk_label_new("");
gtk_widget_set_name(label_x2_eqsecond, "result_label");
gtk_box_pack_start(GTK_BOX(eq_second_box), label_x2_eqsecond, FALSE, FALSE, 0);

gpointer data_eq[5] = {entries_eq2, label_result_eqsecond, label_delta_eqsecond, label_x1_eqsecond, label_x2_eqsecond};

GtkWidget *show_eq_second = gtk_button_new_with_label("Show equation");
gtk_widget_set_name(show_eq_second, "calcmorebutton");
g_signal_connect(show_eq_second, "clicked", G_CALLBACK(on_show_eq_second_button_clicked), data_eq);
gtk_box_pack_start(GTK_BOX(eq_second_box), show_eq_second, FALSE, FALSE, 0);

GtkWidget *eq_second = gtk_button_new_with_label("Square function");
gtk_widget_set_name(eq_second, "morebutton");
g_signal_connect(eq_second, "clicked", G_CALLBACK(on_eq_second_button_clicked), eq_second_box);
gtk_box_pack_start(GTK_BOX(morebox), eq_second, FALSE, FALSE, 0);

///END SQUARE EQUATION

//CUBIC EQUATION
GtkWidget *eq_cubic_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
gtk_widget_set_name(eq_cubic_box, "moreboxes");
gtk_container_add(GTK_CONTAINER(morebox), eq_cubic_box);

GtkWidget *label_eq_cubic = gtk_label_new("Cubic function");
gtk_widget_set_name(label_eq_cubic, "more_label_title");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), label_eq_cubic, FALSE, FALSE, 0);

GtkWidget *label_eq_cubic_formula = gtk_label_new("ax³ + bx² + cx + d = 0");
gtk_widget_set_name(label_eq_cubic_formula, "more_label_title");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), label_eq_cubic_formula, FALSE, FALSE, 0);


GtkWidget *label_a_eqcubic = gtk_label_new("Enter the value for a:");
gtk_widget_set_name(label_a_eqcubic, "enter_label");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), label_a_eqcubic, FALSE, FALSE, 0);

GtkWidget *entrya_eqcubic = gtk_entry_new();
gtk_widget_set_name(entrya_eqcubic, "more_entry");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), entrya_eqcubic, FALSE, FALSE, 0);

GtkWidget *label_b_eqcubic = gtk_label_new("Enter the value for b:");
gtk_widget_set_name(label_b_eqcubic, "enter_label");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), label_b_eqcubic, FALSE, FALSE, 0);

GtkWidget *entryb_eqcubic = gtk_entry_new();
gtk_widget_set_name(entryb_eqcubic, "more_entry");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), entryb_eqcubic, FALSE, FALSE, 0);

GtkWidget *label_c_eqcubic = gtk_label_new("Enter the value for c:");
gtk_widget_set_name(label_c_eqcubic, "enter_label");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), label_c_eqcubic, FALSE, FALSE, 0);

GtkWidget *entryc_eqcubic = gtk_entry_new();
gtk_widget_set_name(entryc_eqcubic, "more_entry");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), entryc_eqcubic, FALSE, FALSE, 0);

GtkWidget *label_d_eqcubic = gtk_label_new("Enter the value for d:");
gtk_widget_set_name(label_d_eqcubic, "enter_label");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), label_d_eqcubic, FALSE, FALSE, 0);

GtkWidget *entryd_eqcubic = gtk_entry_new();
gtk_widget_set_name(entryd_eqcubic, "more_entry");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), entryd_eqcubic, FALSE, FALSE, 0);

// Store the entries in an array
GtkWidget *entries_eq3[4] = {entrya_eqcubic, entryb_eqcubic, entryc_eqcubic, entryd_eqcubic};

GtkWidget *label_result_eqcubic = gtk_label_new("");
gtk_widget_set_name(label_result_eqcubic, "result_label");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), label_result_eqcubic, FALSE, FALSE, 0);

GtkWidget *label_x_eqcubic = gtk_label_new("");
gtk_widget_set_name(label_x_eqcubic, "result_label");
gtk_box_pack_start(GTK_BOX(eq_cubic_box), label_x_eqcubic, FALSE, FALSE, 0);

gpointer data_eqcubic[3] = {entries_eq3, label_result_eqcubic, label_x_eqcubic};

GtkWidget *show_eq_cubic = gtk_button_new_with_label("Show equation");
gtk_widget_set_name(show_eq_cubic, "calcmorebutton");
g_signal_connect(show_eq_cubic, "clicked", G_CALLBACK(on_show_eq_cubic_button_clicked), data_eqcubic);
gtk_box_pack_start(GTK_BOX(eq_cubic_box), show_eq_cubic, FALSE, FALSE, 0);

GtkWidget *eq_cubic = gtk_button_new_with_label("Cubic function");
gtk_widget_set_name(eq_cubic, "morebutton");
g_signal_connect(eq_cubic, "clicked", G_CALLBACK(on_eq_cubic_button_clicked), eq_cubic_box);
gtk_box_pack_start(GTK_BOX(morebox), eq_cubic, FALSE, FALSE, 0);

//END CUBIC EQUATION
    
//ARCHIMEDEAN SPIRAL
GtkWidget *eq_archim_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
gtk_widget_set_name(eq_archim_box, "moreboxes");
gtk_container_add(GTK_CONTAINER(morebox), eq_archim_box);

GtkWidget *label_archim = gtk_label_new("Archimedean spiral");
gtk_widget_set_name(label_archim, "more_label_title");
gtk_box_pack_start(GTK_BOX(eq_archim_box), label_archim, FALSE, FALSE, 0);

GtkWidget *label_archim_formula = gtk_label_new("r = kΘ");
gtk_widget_set_name(label_archim_formula, "more_label_title");
gtk_box_pack_start(GTK_BOX(eq_archim_box), label_archim_formula, FALSE, FALSE, 0);


GtkWidget *label_t_eqarchim = gtk_label_new("Enter the value for t: \ntime");
gtk_widget_set_name(label_t_eqarchim, "enter_label");
gtk_box_pack_start(GTK_BOX(eq_archim_box), label_t_eqarchim, FALSE, FALSE, 0);

GtkWidget *entryt_eqarchim = gtk_entry_new();
gtk_widget_set_name(entryt_eqarchim, "more_entry");
gtk_box_pack_start(GTK_BOX(eq_archim_box), entryt_eqarchim, FALSE, FALSE, 0);

GtkWidget *label_v_eqarchim = gtk_label_new("Enter the value for v: \nvelocity");
gtk_widget_set_name(label_v_eqarchim, "enter_label");
gtk_box_pack_start(GTK_BOX(eq_archim_box), label_v_eqarchim, FALSE, FALSE, 0);

GtkWidget *entryv_eqarchim = gtk_entry_new();
gtk_widget_set_name(entryv_eqarchim, "more_entry");
gtk_box_pack_start(GTK_BOX(eq_archim_box), entryv_eqarchim, FALSE, FALSE, 0);

GtkWidget *label_w_eqarchim = gtk_label_new("Enter the value for w: \nangular velocity");
gtk_widget_set_name(label_w_eqarchim, "enter_label");
gtk_box_pack_start(GTK_BOX(eq_archim_box), label_w_eqarchim, FALSE, FALSE, 0);

GtkWidget *entryw_eqarchim = gtk_entry_new();
gtk_widget_set_name(entryw_eqarchim, "more_entry");
gtk_box_pack_start(GTK_BOX(eq_archim_box), entryw_eqarchim, FALSE, FALSE, 0);

GtkWidget *entries_eqarchim[3] = {entryt_eqarchim, entryv_eqarchim, entryw_eqarchim};

gpointer data_eqarchim[3] = {entries_eqarchim};

GtkWidget *show_archim = gtk_button_new_with_label("Show spiral");
gtk_widget_set_name(show_archim, "calcmorebutton");
g_signal_connect(show_archim, "clicked", G_CALLBACK(on_show_archim_button_clicked), data_eqarchim);
gtk_box_pack_start(GTK_BOX(eq_archim_box), show_archim, FALSE, FALSE, 0);

GtkWidget *eq_archim = gtk_button_new_with_label("Archimedean spiral");
gtk_widget_set_name(eq_archim, "morebutton");
g_signal_connect(eq_archim, "clicked", G_CALLBACK(on_eq_archim_button_clicked), eq_archim_box);
gtk_box_pack_start(GTK_BOX(morebox), eq_archim, FALSE, FALSE, 0);

//END ARCHIMEDEAN SPIRAL   

//EXPONENTIAL FUNCTON
GtkWidget *expo_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
gtk_widget_set_name(expo_box, "moreboxes");
gtk_container_add(GTK_CONTAINER(morebox), expo_box);

GtkWidget *label_expo = gtk_label_new("Exponential function");
gtk_widget_set_name(label_expo, "more_label_title");
gtk_box_pack_start(GTK_BOX(expo_box), label_expo, FALSE, FALSE, 0);

GtkWidget *label_expo_formula = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label_expo_formula), "f(x) = b<sup>X</sup>");
gtk_widget_set_name(label_expo_formula, "more_label_title");
gtk_container_add(GTK_CONTAINER(expo_box), label_expo_formula);

GtkWidget *label2_expo_formula = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label2_expo_formula), "b > 0");
gtk_widget_set_name(label2_expo_formula, "enter_label_expo");
gtk_container_add(GTK_CONTAINER(expo_box), label2_expo_formula);

GtkWidget *label_expo_info = gtk_label_new("Enter up to 6 values for x");
gtk_box_pack_start(GTK_BOX(expo_box), label_expo_info, FALSE, FALSE, 0);

    GtkWidget *xvaluegrid_expo = gtk_grid_new();
    //gtk_widget_set_size_request(xvaluegrid_expo, 200, -1);
    gtk_box_pack_start(GTK_BOX(expo_box), xvaluegrid_expo, TRUE, TRUE, 0);

GtkWidget *label_x1_expo = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label_x1_expo ), "X <sub>1</sub>");
gtk_widget_set_name(label_x1_expo , "enter_label_expo");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), label_x1_expo, 0, 0, 1, 1);

GtkWidget *entryx1_expo = gtk_entry_new();
gtk_widget_set_name(entryx1_expo, "more_entry");
gtk_widget_set_size_request(entryx1_expo, 50, -1);
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), entryx1_expo, 1, 0, 1, 1);

GtkWidget *label_x2_expo = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label_x2_expo ), "X <sub>2</sub>");
gtk_widget_set_name(label_x2_expo , "enter_label_expo");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), label_x2_expo, 0, 1, 1, 1);

GtkWidget *entryx2_expo = gtk_entry_new();
gtk_widget_set_name(entryx2_expo, "more_entry");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), entryx2_expo, 1, 1, 1, 1);

GtkWidget *label_x3_expo = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label_x3_expo ), "X <sub>3</sub>");
gtk_widget_set_name(label_x3_expo , "enter_label_expo");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), label_x3_expo, 0, 2, 1, 1);

GtkWidget *entryx3_expo = gtk_entry_new();
gtk_widget_set_name(entryx3_expo, "more_entry");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), entryx3_expo, 1, 2, 1, 1);

GtkWidget *label_x4_expo = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label_x4_expo ), "X <sub>4</sub>");
gtk_widget_set_name(label_x4_expo , "enter_label_expo");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), label_x4_expo, 0, 3, 1, 1);

GtkWidget *entryx4_expo = gtk_entry_new();
gtk_widget_set_name(entryx4_expo, "more_entry");
gtk_widget_set_size_request(entryx4_expo, 50, -1);
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), entryx4_expo, 1, 3, 1, 1);

GtkWidget *label_x5_expo = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label_x5_expo ), "X <sub>5</sub>");
gtk_widget_set_name(label_x5_expo , "enter_label_expo");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), label_x5_expo, 0, 4, 1, 1);

GtkWidget *entryx5_expo = gtk_entry_new();
gtk_widget_set_name(entryx5_expo, "more_entry");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), entryx5_expo, 1, 4, 1, 1);

GtkWidget *label_x6_expo = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label_x6_expo ), "X <sub>6</sub>");
gtk_widget_set_name(label_x6_expo , "enter_label_expo");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), label_x6_expo, 0, 5, 1, 1);

GtkWidget *entryx6_expo = gtk_entry_new();
gtk_widget_set_name(entryx6_expo, "more_entry");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), entryx6_expo, 1, 5, 1, 1);

GtkWidget *label_b_expo = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label_b_expo ), "b");
gtk_widget_set_name(label_b_expo , "enter_label_expo");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), label_b_expo, 0, 6, 1, 1);

GtkWidget *entryb_expo = gtk_entry_new();
gtk_widget_set_name(entryb_expo, "more_entry");
gtk_grid_attach(GTK_GRID(xvaluegrid_expo), entryb_expo, 1, 6, 1, 1);


GtkWidget *entries_expo[7] = {entryx1_expo, entryx2_expo, entryx3_expo, entryx4_expo, entryx5_expo, entryx6_expo, entryb_expo,};

gpointer data_expo[1] = {entries_expo};

GtkWidget *show_expo = gtk_button_new_with_label("Show result");
gtk_widget_set_name(show_expo, "calcmorebutton");
g_signal_connect(show_expo, "clicked", G_CALLBACK(on_show_expo_button_clicked), data_expo);
gtk_box_pack_start(GTK_BOX(expo_box), show_expo, FALSE, FALSE, 0);

GtkWidget *eq_expo = gtk_button_new_with_label("Exponential function");
gtk_widget_set_name(eq_expo, "morebutton");
g_signal_connect(eq_expo, "clicked", G_CALLBACK(on_eq_expo_button_clicked), expo_box);
gtk_box_pack_start(GTK_BOX(morebox), eq_expo, FALSE, FALSE, 0);

//END EXPONENTIAL FUNCTION

//START PROPORTIONS
GtkWidget *proportions_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
gtk_widget_set_name(proportions_box, "moreboxes");
gtk_container_add(GTK_CONTAINER(morebox), proportions_box);

GtkWidget *label_proportions = gtk_label_new("Proportions");
gtk_widget_set_name(label_proportions, "more_label_title");
gtk_box_pack_start(GTK_BOX(proportions_box), label_proportions, FALSE, FALSE, 0);

GtkWidget *label_proportions_formula = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label_proportions_formula), "a : b = c : d");
gtk_widget_set_name(label_proportions_formula, "more_label_title");
gtk_container_add(GTK_CONTAINER(proportions_box), label_proportions_formula);

GtkWidget *label_proportions_info = gtk_label_new("Enter x for one value");
gtk_box_pack_start(GTK_BOX(proportions_box), label_proportions_info, FALSE, FALSE, 0);

GtkWidget *label_a_proportions = gtk_label_new("Enter the value for a:");
gtk_widget_set_name(label_a_proportions, "enter_label");
gtk_box_pack_start(GTK_BOX(proportions_box), label_a_proportions, FALSE, FALSE, 0);

GtkWidget *entrya_proportions = gtk_entry_new();
gtk_widget_set_name(entrya_proportions, "more_entry");
gtk_box_pack_start(GTK_BOX(proportions_box), entrya_proportions, FALSE, FALSE, 0);

GtkWidget *label_b_proportions = gtk_label_new("Enter the value for b:");
gtk_widget_set_name(label_b_proportions, "enter_label");
gtk_box_pack_start(GTK_BOX(proportions_box), label_b_proportions, FALSE, FALSE, 0);

GtkWidget *entryb_proportions = gtk_entry_new();
gtk_widget_set_name(entryb_proportions, "more_entry");
gtk_box_pack_start(GTK_BOX(proportions_box), entryb_proportions, FALSE, FALSE, 0);

GtkWidget *label_c_proportions = gtk_label_new("Enter the value for c:");
gtk_widget_set_name(label_c_proportions, "enter_label");
gtk_box_pack_start(GTK_BOX(proportions_box), label_c_proportions, FALSE, FALSE, 0);

GtkWidget *entryc_proportions = gtk_entry_new();
gtk_widget_set_name(entryc_proportions, "more_entry");
gtk_box_pack_start(GTK_BOX(proportions_box), entryc_proportions, FALSE, FALSE, 0);

GtkWidget *label_d_proportions = gtk_label_new("Enter the value for d:");
gtk_widget_set_name(label_d_proportions, "enter_label");
gtk_box_pack_start(GTK_BOX(proportions_box), label_d_proportions, FALSE, FALSE, 0);

GtkWidget *entryd_proportions = gtk_entry_new();
gtk_widget_set_name(entryd_proportions, "more_entry");
gtk_box_pack_start(GTK_BOX(proportions_box), entryd_proportions, FALSE, FALSE, 0);

// RESULT
label_result_proportions = gtk_label_new("");
gtk_widget_set_name(label_result_proportions, "result_label");
gtk_box_pack_start(GTK_BOX(proportions_box), label_result_proportions, FALSE, FALSE, 0);

// X
label_x_proportions = gtk_label_new("");
gtk_widget_set_name(label_x_proportions, "result_label");
gtk_box_pack_start(GTK_BOX(proportions_box), label_x_proportions, FALSE, FALSE, 0);

GtkWidget *entries_proportions[4] = {entrya_proportions, entryb_proportions, entryc_proportions, entryd_proportions};

gpointer data_proportions[3] = {entries_proportions, label_x_proportions, label_result_proportions};

GtkWidget *show_proportions = gtk_button_new_with_label("Show result");
gtk_widget_set_name(show_proportions, "calcmorebutton");
g_signal_connect(show_proportions, "clicked", G_CALLBACK(on_show_proportions_button_clicked), data_proportions);
gtk_box_pack_start(GTK_BOX(proportions_box), show_proportions, FALSE, FALSE, 0);

GtkWidget *eq_proportions = gtk_button_new_with_label("Proportions");
gtk_widget_set_name(eq_proportions, "morebutton");
g_signal_connect(eq_proportions, "clicked", G_CALLBACK(on_eq_proportions_button_clicked), proportions_box);
gtk_box_pack_start(GTK_BOX(morebox), eq_proportions, FALSE, FALSE, 0);
//END PROPORTIONS

//LINEAR INEQUALITY
GtkWidget *ineq2_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
gtk_widget_set_name(ineq2_box, "moreboxes");
gtk_container_add(GTK_CONTAINER(morebox), ineq2_box);

GtkWidget *label_ineq2 = gtk_label_new("Inequality x y");
gtk_widget_set_name(label_ineq2, "more_label_title");
gtk_box_pack_start(GTK_BOX(ineq2_box), label_ineq2, FALSE, FALSE, 0);

GtkWidget *label_ineq2_formula = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label_ineq2_formula), "ax + by \u2264 c ");
gtk_widget_set_name(label_ineq2_formula, "more_label_title");
gtk_container_add(GTK_CONTAINER(ineq2_box), label_ineq2_formula);

GtkWidget *label_ineq2_formula2 = gtk_label_new(nullptr);
gtk_label_set_markup(GTK_LABEL(label_ineq2_formula2), "ax + by \u2265 c ");
gtk_widget_set_name(label_ineq2_formula2, "more_label_title");
gtk_container_add(GTK_CONTAINER(ineq2_box), label_ineq2_formula2);

GtkWidget *label_a_ineq2 = gtk_label_new("Enter the value for a:");
gtk_widget_set_name(label_a_ineq2, "enter_label");
gtk_box_pack_start(GTK_BOX(ineq2_box), label_a_ineq2, FALSE, FALSE, 0);
//ENTRIES
GtkWidget *entrya_ineq2 = gtk_entry_new();
gtk_widget_set_name(entrya_ineq2, "more_entry");
gtk_box_pack_start(GTK_BOX(ineq2_box), entrya_ineq2, FALSE, FALSE, 0);

GtkWidget *label_b_ineq2 = gtk_label_new("Enter the value for b:");
gtk_widget_set_name(label_b_ineq2, "enter_label");
gtk_box_pack_start(GTK_BOX(ineq2_box), label_b_ineq2, FALSE, FALSE, 0);

GtkWidget *entryb_ineq2 = gtk_entry_new();
gtk_widget_set_name(entryb_ineq2, "more_entry");
gtk_box_pack_start(GTK_BOX(ineq2_box), entryb_ineq2, FALSE, FALSE, 0);

GtkWidget *label_c_ineq2 = gtk_label_new("Enter the value for c:");
gtk_widget_set_name(label_c_ineq2, "enter_label");
gtk_box_pack_start(GTK_BOX(ineq2_box), label_c_ineq2, FALSE, FALSE, 0);

GtkWidget *entryc_ineq2 = gtk_entry_new();
gtk_widget_set_name(entryc_ineq2, "more_entry");
gtk_box_pack_start(GTK_BOX(ineq2_box), entryc_ineq2, FALSE, FALSE, 0);

// RESULT
label_result_ineq2 = gtk_label_new("");
gtk_widget_set_name(label_result_ineq2, "result_label");
gtk_box_pack_start(GTK_BOX(ineq2_box), label_result_ineq2, FALSE, FALSE, 0);

// X
label_x_ineq2 = gtk_label_new("");
gtk_widget_set_name(label_x_ineq2, "result_label");
gtk_box_pack_start(GTK_BOX(ineq2_box), label_x_ineq2, FALSE, FALSE, 0);
// Y
label_y_ineq2 = gtk_label_new("");
gtk_widget_set_name(label_y_ineq2, "result_label");
gtk_box_pack_start(GTK_BOX(ineq2_box), label_y_ineq2, FALSE, FALSE, 0);
// X2
label_x2_ineq2 = gtk_label_new("");
gtk_widget_set_name(label_x2_ineq2, "result_label");
gtk_box_pack_start(GTK_BOX(ineq2_box), label_x2_ineq2, FALSE, FALSE, 0);
// Y2
label_y2_ineq2 = gtk_label_new("");
gtk_widget_set_name(label_y2_ineq2, "result_label");
gtk_box_pack_start(GTK_BOX(ineq2_box), label_y2_ineq2, FALSE, FALSE, 0);

GtkWidget *entries_ineq2[3] = {entrya_ineq2, entryb_ineq2, entryc_ineq2};

gpointer data_ineq2[6] = {entries_ineq2, label_result_ineq2, label_x_ineq2, label_y_ineq2, label_x2_ineq2, label_y2_ineq2};

GtkWidget *show_ineq2 = gtk_button_new_with_label("Show result for \u2264");
gtk_widget_set_name(show_ineq2, "calcmorebutton");
g_signal_connect(show_ineq2, "clicked", G_CALLBACK(on_show_ineq2_button_clicked), data_ineq2);
gtk_box_pack_start(GTK_BOX(ineq2_box), show_ineq2, FALSE, FALSE, 0);

GtkWidget *eq_ineq2 = gtk_button_new_with_label("Inequality x y");
gtk_widget_set_name(eq_ineq2, "morebutton");
g_signal_connect(eq_ineq2, "clicked", G_CALLBACK(on_eq_ineq2_button_clicked), ineq2_box);
gtk_box_pack_start(GTK_BOX(morebox), eq_ineq2, FALSE, FALSE, 0);
//END LINEAR INEQUALITY

//SIGNAL
g_signal_connect(equal, "clicked", G_CALLBACK(on_equal_button_clicked), entry_expression);
g_object_set_data(G_OBJECT(equal), "entry_result", entry_result); // Pass entry_result to the equal function

// Create the box
GtkWidget *calc_notebox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_widget_set_size_request(calc_notebox, 500, 1);
gtk_box_pack_end(GTK_BOX(calculatorbar), calc_notebox, TRUE, TRUE, 0);
gtk_widget_set_hexpand(calc_notebox, TRUE);
gtk_widget_set_vexpand(calc_notebox, TRUE);

// Create a scrolled window for better text handling
GtkWidget *calc_notebox_scrolled_window = gtk_scrolled_window_new(NULL, NULL);
gtk_widget_set_hexpand(calc_notebox_scrolled_window, TRUE);
gtk_widget_set_vexpand(calc_notebox_scrolled_window, TRUE);
gtk_container_add(GTK_CONTAINER(calc_notebox), calc_notebox_scrolled_window);

// Create a text view
GtkWidget *calc_notebox_text_view = gtk_text_view_new();
gtk_widget_set_name(calc_notebox_text_view, "calc_notebox_text_view");
gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(calc_notebox_text_view), GTK_WRAP_WORD); // Word wrap for better formatting
gtk_container_add(GTK_CONTAINER(calc_notebox_scrolled_window), calc_notebox_text_view);
gtk_text_view_set_editable(GTK_TEXT_VIEW(calc_notebox_text_view), TRUE);
gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(calc_notebox_text_view), TRUE);



//END CALCULATOR



//COLOR PICKER
GtkWidget *colorpickerbar = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_widget_set_name(colorpickerbar, "colorpickerbar");
gtk_widget_set_size_request(colorpickerbar, 300, 1); 
gtk_box_pack_start(GTK_BOX(hbox), colorpickerbar, FALSE, TRUE, 0);

GtkWidget *color_label = gtk_label_new("COLOR PICKER");
gtk_widget_set_halign(color_label, GTK_ALIGN_CENTER);
gtk_widget_set_name(color_label, "web_history_label");
gtk_box_pack_start(GTK_BOX(colorpickerbar), color_label, false, false, 0);

GtkWidget *rgb_entry = gtk_entry_new();
gtk_widget_set_name(rgb_entry , "colorentry");
gtk_box_pack_start(GTK_BOX(colorpickerbar), rgb_entry, FALSE, TRUE, 0);

GtkWidget *hex_entry = gtk_entry_new();
gtk_widget_set_name(hex_entry, "colorhexentry");
gtk_box_pack_start(GTK_BOX(colorpickerbar), hex_entry, FALSE, TRUE, 0);

    GtkWidget *color_area = gtk_drawing_area_new();
    gtk_widget_set_name(color_area, "colorarea");
    gtk_widget_set_size_request(color_area, 300, 300);
    gtk_box_pack_start(GTK_BOX(colorpickerbar), color_area, FALSE, TRUE, 0);

    GtkWidget *brightness_slider = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, 0.0, 2.0, 0.01);
    gtk_range_set_value(GTK_RANGE(brightness_slider), brightness);
    gtk_box_pack_start(GTK_BOX(colorpickerbar), brightness_slider, FALSE, TRUE, 0);
    g_signal_connect(brightness_slider, "value-changed", G_CALLBACK(on_brightness_changed), color_area);

GtkWidget *space = gtk_entry_new();
gtk_widget_set_name(space, "space");
gtk_box_pack_start(GTK_BOX(colorpickerbar), space, FALSE, TRUE, 0);

// Make the entry not clickable
gtk_widget_set_sensitive(space, FALSE);


    GtkWidget *colorbw_area = gtk_drawing_area_new();
    gtk_widget_set_name(colorbw_area, "colorbw_area");
    gtk_widget_set_size_request(colorbw_area, 300, 25);
    gtk_box_pack_start(GTK_BOX(colorpickerbar), colorbw_area, FALSE, TRUE, 0);
    
    gtk_widget_add_events(color_area, GDK_BUTTON_PRESS_MASK);
    g_signal_connect(color_area, "draw", G_CALLBACK(on_draw), NULL);
    g_signal_connect(color_area, "button-press-event", G_CALLBACK(on_color_selected), NULL);

    gtk_widget_add_events(colorbw_area, GDK_BUTTON_PRESS_MASK);
    g_signal_connect(colorbw_area, "draw", G_CALLBACK(on_drawbw), NULL);
    g_signal_connect(colorbw_area, "button-press-event", G_CALLBACK(on_colorbw_selected), NULL);

    g_object_set_data(G_OBJECT(color_area), "rgb_entry", rgb_entry);
    g_object_set_data(G_OBJECT(color_area), "hex_entry", hex_entry);
    g_object_set_data(G_OBJECT(colorbw_area), "rgb_entry", rgb_entry);
    g_object_set_data(G_OBJECT(colorbw_area), "hex_entry", hex_entry);

GtkWidget *space2 = gtk_entry_new();
gtk_widget_set_name(space2, "space");
gtk_box_pack_start(GTK_BOX(colorpickerbar), space2, FALSE, TRUE, 0);
    
    
    GtkWidget *showcolor = gtk_drawing_area_new();
    gtk_widget_set_name(showcolor, "showcolor");
    gtk_widget_set_size_request(showcolor, 200, 150);
    gtk_box_pack_start(GTK_BOX(colorpickerbar), showcolor, FALSE, TRUE, 0);    
    g_signal_connect(G_OBJECT(showcolor), "draw", G_CALLBACK(drawcolor_callback), NULL);
    g_object_set_data(G_OBJECT(color_area), "showcolor", showcolor);
    g_object_set_data(G_OBJECT(colorbw_area), "showcolor", showcolor);
//END COLOR PICKER    
  
//SPEED TEST BAR
GtkWidget *speedtestbar = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_widget_set_name(speedtestbar, "speedtestbar");
gtk_widget_set_size_request(colorpickerbar, 300, 1); 
gtk_box_pack_start(GTK_BOX(hbox), speedtestbar, FALSE, TRUE, 0);

GtkWidget* svbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
gtk_container_add(GTK_CONTAINER(speedtestbar), svbox);

GtkWidget* speedimg = gtk_button_new_with_label("");
gtk_widget_set_name(speedimg, "speedimg");
gtk_widget_set_halign(speedimg, GTK_ALIGN_CENTER);
gtk_box_pack_start(GTK_BOX(svbox), speedimg, FALSE, FALSE, 0);


// Download Entry
GtkWidget *download_title = gtk_label_new("Download Speed");
gtk_widget_set_name(download_title, "inside_menu_title");
gtk_box_pack_start(GTK_BOX(svbox), download_title, FALSE, FALSE, 5);

GtkWidget* download_entry = gtk_entry_new();
gtk_entry_set_placeholder_text(GTK_ENTRY(download_entry), "");
gtk_widget_set_name(download_entry, "colorentry");
gtk_editable_set_editable(GTK_EDITABLE(download_entry), FALSE);
gtk_box_pack_start(GTK_BOX(svbox), download_entry, FALSE, FALSE, 5);

// Upload Entry
GtkWidget *upload_title = gtk_label_new("Upload Speed");
gtk_widget_set_name(upload_title, "inside_menu_title");
gtk_box_pack_start(GTK_BOX(svbox), upload_title, FALSE, FALSE, 5);

GtkWidget* upload_entry = gtk_entry_new();
gtk_entry_set_placeholder_text(GTK_ENTRY(upload_entry), "");
gtk_widget_set_name(upload_entry, "colorentry");
gtk_editable_set_editable(GTK_EDITABLE(upload_entry), FALSE);
gtk_box_pack_start(GTK_BOX(svbox), upload_entry, FALSE, FALSE, 5);

// Spinner for Loading Indicator
GtkWidget* spinner = gtk_spinner_new();
gtk_widget_set_name(spinner, "spinner");
gtk_box_pack_start(GTK_BOX(svbox), spinner, FALSE, FALSE, 5);

// Speed Test Button
GtkWidget* speedbutton = gtk_button_new_with_label("Start Test Speed");
gtk_widget_set_name(speedbutton, "add_button");
gtk_box_pack_start(GTK_BOX(svbox), speedbutton, FALSE, FALSE, 0);

// Combine into one structure for callback
struct SpeedTestData {
    GtkEntry* download_entry;
    GtkEntry* upload_entry;
    GtkWidget* spinner;
};

SpeedTestData* data_speed = new SpeedTestData{
    GTK_ENTRY(download_entry),
    GTK_ENTRY(upload_entry),
    spinner,
};

// Signal Connection for Button Click
g_signal_connect(speedbutton, "clicked", G_CALLBACK(on_test_speed_clicked), data_speed);

//END SPEED TEST
    
//
  
//INSTRUMENTSBAR 
    GtkWidget *instrumentbar = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_widget_set_name(instrumentbar, "instrumentbar");
    gtk_widget_set_size_request(instrumentbar, 48, 1); 
    gtk_box_pack_start(GTK_BOX(hbox), instrumentbar, false, true, 0);

    GtkWidget *instrumentbargrid = gtk_grid_new();
    gtk_container_set_border_width(GTK_CONTAINER(instrumentbargrid), 5);
    gtk_box_pack_start(GTK_BOX(instrumentbar), instrumentbargrid, TRUE, TRUE, 0);
    
    // NOTE button
    GtkWidget *note_button = gtk_button_new_with_label("");
    gtk_widget_set_name(note_button, "note_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(note_button), "Notebook"); 
    gtk_grid_attach(GTK_GRID(instrumentbargrid), note_button, 0, 0, 1, 1); 
    g_signal_connect(note_button, "clicked", G_CALLBACK(on_note_button_clicked), notebar);    
    // COLOR PICKER button
    GtkWidget *colorpicker_button = gtk_button_new_with_label("");
    gtk_widget_set_name(colorpicker_button, "colorpicker_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(colorpicker_button), "Color picker"); 
    gtk_grid_attach(GTK_GRID(instrumentbargrid), colorpicker_button, 0, 1, 1, 1); 
    g_signal_connect(colorpicker_button, "clicked", G_CALLBACK(on_colorpicker_button_clicked), colorpickerbar);   
    // CALCULATOR button
    GtkWidget *calculator_button = gtk_button_new_with_label("");
    gtk_widget_set_name(calculator_button, "calculator_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(calculator_button), "Calculator"); 
    gtk_grid_attach(GTK_GRID(instrumentbargrid), calculator_button, 0, 2, 1, 1);
    g_signal_connect(calculator_button, "clicked", G_CALLBACK(on_calculator_button_clicked), calculatorbar);    
    //CONSOLE button
    GtkWidget *console_button = gtk_button_new_with_label("");
    gtk_widget_set_name(console_button, "console_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(console_button), "Console"); 
    gtk_grid_attach(GTK_GRID(instrumentbargrid), console_button, 0, 3, 1, 1); 
    g_signal_connect(console_button, "clicked", G_CALLBACK(console_button_clicked), notebook);
    //SPEED TEST button
    GtkWidget *speedtest_button = gtk_button_new_with_label("");
    gtk_widget_set_name(speedtest_button, "speedtest_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(speedtest_button), "Speed Test"); 
    gtk_grid_attach(GTK_GRID(instrumentbargrid), speedtest_button, 0, 4, 1, 1); 
    g_signal_connect(speedtest_button, "clicked", G_CALLBACK(on_speedtest_button_clicked), speedtestbar);

//END INSTRUMENT BAR  



// TOOLBAR
    GtkWidget *toolbar = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_widget_set_size_request(toolbar, 48, -1);     
    gtk_box_pack_start(GTK_BOX(hbox), toolbar, false, true, 0);  

    GtkWidget *grid = gtk_grid_new();
    gtk_container_set_border_width(GTK_CONTAINER(grid), 13);
    gtk_box_pack_start(GTK_BOX(toolbar), grid, TRUE, TRUE, 0);
   
    // MENU button
    GtkWidget *menu_button = gtk_button_new_with_label("");
    gtk_widget_set_name(menu_button, "menu_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(menu_button), "Menu"); 
    gtk_grid_attach(GTK_GRID(grid), menu_button, 0, 0, 1, 1); // Column 0, Row 1
    g_signal_connect(menu_button, "clicked", G_CALLBACK(on_menu_button_clicked), menubar);   
    // PLUS BUTTON
    GtkWidget *new_button = gtk_button_new_with_label("+");
    gtk_widget_set_name(new_button, "button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(new_button), "Add tab"); 
    gtk_grid_attach(GTK_GRID(grid), new_button, 0, 1, 1, 1);
    g_signal_connect(new_button, "clicked", G_CALLBACK(on_new_button_clicked), notebook);   
    // CLOSETAB BUTTON
    GtkWidget *closetab_button = gtk_button_new_with_label("-");
    gtk_widget_set_name(closetab_button, "closetab_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(closetab_button), "Close tab"); 
    gtk_grid_attach(GTK_GRID(grid), closetab_button, 0, 2, 1, 1);   
    g_signal_connect(closetab_button, "clicked", G_CALLBACK(close_active_tab), notebook); 
    // RELOAD button
    GtkWidget *reload_button = gtk_button_new_with_label("");
    gtk_widget_set_name(reload_button, "reload_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(reload_button), "Reload back"); 
    gtk_grid_attach(GTK_GRID(grid), reload_button, 0, 3, 1, 1);
    g_signal_connect(reload_button, "clicked", G_CALLBACK(on_reload_button_clicked), notebook);   
    // BACK button
    GtkWidget *back_button = gtk_button_new_with_label("<");
    gtk_widget_set_name(back_button, "back_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(back_button), "Go back");              
    gtk_grid_attach(GTK_GRID(grid), back_button, 0, 4, 1, 1); 
    g_signal_connect(back_button, "clicked", G_CALLBACK(on_back_button_clicked), notebook);   
    // FORWARD button
    GtkWidget *forward_button = gtk_button_new_with_label(">");
    gtk_widget_set_name(forward_button, "forward_button"); 
    gtk_widget_set_tooltip_text(GTK_WIDGET(forward_button), "Go forward");              
    gtk_grid_attach(GTK_GRID(grid), forward_button, 0, 5, 1, 1); 
    g_signal_connect(forward_button, "clicked", G_CALLBACK(on_forward_button_clicked), notebook);
    // SEARCH button
    GtkWidget *search_button = gtk_button_new_with_label("");
    gtk_widget_set_name(search_button, "search_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(search_button), "Search Bar");          
    gtk_grid_attach(GTK_GRID(grid), search_button, 0, 6, 1, 1); 
    g_signal_connect(search_button, "clicked", G_CALLBACK(on_search_button_clicked), entry);    
    // SEARCH ENGINE button
    GtkWidget *searchengine_button = gtk_button_new_with_label("");
    gtk_widget_set_name(searchengine_button, "searchengine_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(searchengine_button), "Search Engines");          
    gtk_grid_attach(GTK_GRID(grid), searchengine_button, 0, 7, 1, 1); 
    g_signal_connect(searchengine_button, "clicked", G_CALLBACK(on_searchengine_button_clicked), searchenginebar);
    // TOOLS button
    GtkWidget *tools_button = gtk_button_new_with_label("");
    gtk_widget_set_name(tools_button, "tools_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(tools_button), "Tools");          
    gtk_grid_attach(GTK_GRID(grid), tools_button, 0, 8, 1, 1); 
    g_signal_connect(tools_button, "clicked", G_CALLBACK(on_tools_button_clicked), instrumentbar);    
    // BOOKMARK button
    GtkWidget *openbookmark_button = gtk_button_new_with_label("");
    gtk_widget_set_name(openbookmark_button, "openbookmark_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(openbookmark_button), "Bookmarks");      
    gtk_grid_attach(GTK_GRID(grid), openbookmark_button, 0, 9, 1, 1); 
    g_signal_connect(openbookmark_button, "clicked", G_CALLBACK(on_openbookmark_button_clicked), bookmarkbar);  
    // DOWNLOAD button
    GtkWidget *download_button = gtk_button_new_with_label("");
    gtk_widget_set_name(download_button, "download_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(download_button), "Download");    
    gtk_grid_attach(GTK_GRID(grid), download_button, 0, 10, 1, 1); 
    g_signal_connect(download_button, "clicked", G_CALLBACK(on_download_button_clicked), downloadbar); 
    // HISTORY button
    GtkWidget *openhistory_button = gtk_button_new_with_label("");
    gtk_widget_set_name(openhistory_button, "openhistory_button");
    gtk_widget_set_tooltip_text(GTK_WIDGET(openhistory_button), "History");
    gtk_grid_attach(GTK_GRID(grid), openhistory_button, 0, 11, 1, 1); 
    g_signal_connect(openhistory_button, "clicked", G_CALLBACK(on_openhistory_button_clicked), historybar);   
    
    gtk_widget_show_all(toolbar);
////END TOOLBAR


    //FUNCTION ACTIVATION (must stay here)
    g_signal_connect(notebook, "switch-page", G_CALLBACK(on_switch_page), entry);
    
   // Set user data to be accessed in callback
    g_object_set_data(G_OBJECT(addbookmark_button), "bookmarkbar", bookmarkbar);
    load_bookmarks(bookmarkbar);


    // SHOW/HIDE
    gtk_widget_show_all(window);

    gtk_widget_hide(menubar); 
    gtk_widget_hide(searchenginebar); 
    gtk_widget_hide(notebar); 
    gtk_widget_hide(bookmarkbar); 
    gtk_widget_hide(instrumentbar);
    gtk_widget_hide(downloadbar);
    gtk_widget_hide(colorpickerbar);
    gtk_widget_hide(speedtestbar);    
    gtk_widget_hide(historybar); 
    gtk_widget_hide(calculatorbar); 

    //gtk_widget_hide(current_download_box);
    //gtk_widget_hide(themebox);   
    //gtk_widget_hide(privacybox);   
    
    
    //calculator
       //gtk_widget_hide(morebox);
    gtk_widget_hide(parabolabox);
    gtk_widget_hide(parabola_x_box);
    gtk_widget_hide(line_box);
    gtk_widget_hide(eq_second_box);
    gtk_widget_hide(eq_cubic_box);
    gtk_widget_hide(eq_archim_box);
    gtk_widget_hide(expo_box);
    gtk_widget_hide(proportions_box);
    gtk_widget_hide(ineq2_box);

    
    gtk_main();

    return 0;
}

 
