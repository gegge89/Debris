#include <gtk/gtk.h>

// Function to load CSS file
static void load_css(const char *css_file) {
    GtkCssProvider *provider = gtk_css_provider_new();
    gboolean success = gtk_css_provider_load_from_path(provider, css_file, NULL);
    if (success) {
        gtk_style_context_add_provider_for_screen(gdk_screen_get_default(),
                                                  GTK_STYLE_PROVIDER(provider),
                                                  GTK_STYLE_PROVIDER_PRIORITY_USER);
        g_print("CSS file '%s' loaded successfully.\n", css_file);
    } else {
        g_warning("Failed to load CSS file '%s'.\n", css_file);
    }
    g_object_unref(provider);
}

static gboolean on_button_press_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
    if (event->button == GDK_BUTTON_PRIMARY) {
        gtk_window_begin_move_drag(GTK_WINDOW(widget), event->button, event->x_root, event->y_root, event->time);
}
    return TRUE;
}

// Close Button Handler
static void on_close_button_clicked(GtkWidget *button, gpointer user_data) {
    GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(user_data),
                                               GTK_DIALOG_MODAL,
                                               GTK_MESSAGE_QUESTION,
                                               GTK_BUTTONS_NONE,
                                               "Are you sure?");
    gtk_widget_set_name(dialog, "dialog");

    // Add custom buttons with correct names
    GtkWidget *yes_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Yes", GTK_RESPONSE_YES);
    gtk_widget_set_name(yes_button, "yes_button");
    GtkWidget *no_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "No", GTK_RESPONSE_NO);
    gtk_widget_set_name(no_button, "no_button");

    // Set default response to "No"
    gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_NO);

    gint response = gtk_dialog_run(GTK_DIALOG(dialog));
    if (response == GTK_RESPONSE_YES) {
        gtk_widget_destroy(GTK_WIDGET(user_data));
    }

    gtk_widget_destroy(dialog);
}

// Minimize Button Handler
static void on_minimize_button_clicked(GtkWidget *widget, gpointer data) {
    gtk_window_iconify(GTK_WINDOW(data));
}

// Maximize Button Handler
static void on_maximize_button_clicked(GtkWidget *widget, gpointer data) {
    gtk_window_maximize(GTK_WINDOW(data));
}

// Search Button Handler
static void on_search_button_clicked(GtkWidget *search_button, gpointer user_data) {
    GtkWidget *entry = GTK_WIDGET(user_data);

    // Toggle visibility of the address bar
    if (gtk_widget_get_visible(entry)) {
        gtk_widget_hide(entry);
    } else {
        gtk_widget_show(entry);
        gtk_widget_grab_focus(entry);
    }
}

int main(int argc, char *argv[]) {
    gtk_init(&argc, &argv);

    // Load the CSS file
    load_css("/usr/share/debris/css/dark_theme.css");

    // Create a new window
    GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window), "BROWS");
    //gtk_window_set_default_size(GTK_WINDOW(window), 500, 600);
    gtk_window_set_decorated(GTK_WINDOW(window), false);
    g_signal_connect(window, "button-press-event", G_CALLBACK(on_button_press_event), NULL);
    g_signal_connect(window, "destroy", G_CALLBACK(gtk_main_quit), NULL);

    // Create a vertical box to hold the address bar, notebook, and button
    GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
    gtk_widget_set_name(vbox, "vbox");
    gtk_container_add(GTK_CONTAINER(window), vbox);

    // Create a box to hold the address bar and main box
    GtkWidget *searchbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_widget_set_name(searchbox, "searchbox");
    gtk_container_add(GTK_CONTAINER(vbox), searchbox);

    // Headbar
    GtkWidget *headbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_widget_set_size_request(headbar, 1, 40);
    gtk_widget_set_name(headbar, "headbar");
    gtk_box_pack_start(GTK_BOX(searchbox), headbar, FALSE, TRUE, 0);

    // Label for Headbar
    GtkWidget *label = gtk_label_new("ABOUT");
    gtk_widget_set_name(label, "debris_label");
    gtk_box_pack_start(GTK_BOX(headbar), label, TRUE, TRUE, 0);

    // Button container
    GtkWidget *button_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_widget_set_name(button_box, "button_box");
    gtk_box_pack_end(GTK_BOX(headbar), button_box, FALSE, TRUE, 0);

    // Close button
    GtkWidget *close_button = gtk_button_new_with_label("");
    gtk_widget_set_name(close_button, "close_button");
    gtk_widget_set_tooltip_text(close_button, "Close");
    gtk_widget_set_size_request(close_button, 24, 24);
    g_signal_connect(close_button, "clicked", G_CALLBACK(on_close_button_clicked), window);
    gtk_box_pack_end(GTK_BOX(button_box), close_button, FALSE, TRUE, 0);

    // Maximize button
    /*
    GtkWidget *maximize_button = gtk_button_new_with_label("");
    gtk_widget_set_name(maximize_button, "maximize_button");
    gtk_widget_set_tooltip_text(maximize_button, "Maximize");
    gtk_widget_set_size_request(maximize_button, 24, 24);
    g_signal_connect(maximize_button, "clicked", G_CALLBACK(on_maximize_button_clicked), window);
    gtk_box_pack_end(GTK_BOX(button_box), maximize_button, FALSE, TRUE, 0);*/

    // Minimize button
    GtkWidget *minimize_button = gtk_button_new_with_label("");
    gtk_widget_set_name(minimize_button, "minimize_button");
    gtk_widget_set_tooltip_text(minimize_button, "Minimize");
    gtk_widget_set_size_request(minimize_button, 24, 24);
    g_signal_connect(minimize_button, "clicked", G_CALLBACK(on_minimize_button_clicked), window);
    gtk_box_pack_end(GTK_BOX(button_box), minimize_button, FALSE, TRUE, 0);

GtkWidget *mainpage_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
gtk_widget_set_name(mainpage_box, "mainpage_box_about");
gtk_widget_set_halign(mainpage_box, GTK_ALIGN_FILL);
gtk_widget_set_valign(mainpage_box, GTK_ALIGN_FILL);
gtk_box_pack_start(GTK_BOX(vbox), mainpage_box, FALSE, FALSE, 0);

// Create and add the mainpage_icon button
GtkWidget *mainpage_icon = gtk_button_new_with_label("");
gtk_widget_set_size_request(mainpage_icon, 250, 250);
gtk_widget_set_name(mainpage_icon, "mainpage_icon");
gtk_widget_set_halign(mainpage_icon, GTK_ALIGN_CENTER);
gtk_widget_set_valign(mainpage_icon, GTK_ALIGN_CENTER);
gtk_box_pack_start(GTK_BOX(mainpage_box), mainpage_icon, FALSE, FALSE, 0);


    GtkWidget *version_button = gtk_button_new_with_label("VERSION 1.0");
    //gtk_widget_set_name(version_button, "version_button");
    gtk_widget_set_name(version_button, "mainmenu_open_button");
    gtk_widget_set_halign(version_button, GTK_ALIGN_CENTER);
    gtk_box_pack_start(GTK_BOX(mainpage_box), version_button, false, false, 0);

    
    
    gtk_widget_show_all(window);
    gtk_main();

    return 0;
}


 
