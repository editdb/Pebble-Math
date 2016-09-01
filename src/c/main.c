#include <pebble.h>
#include "SmallMaths.h"

static Window *main_window;
//static GRect window_bounds;

//-------------------------------------------------------------------------------------------------------- 
void init() {
  // Create main Window element and assign to pointer
  main_window = window_create();

  // Set handlers to manage the elements inside the Window
//  window_set_window_handlers(main_window, (WindowHandlers) {
//    .load = main_window_load,
//    .unload = main_window_unload
//  });

  // Show the Window on the watch, with animated=true
  window_stack_push(main_window, true);
  
  float x = 0.2;
  float y = 0.0;

//  for (int i=0; i<20; i++) {
//    y = sm_ln(x);
//    printf("sm_ln(%i/1000)=%i/1000", (int)(1000.0*x), (int)(1000.0*y));
//    x *= 2;
//  }

  x = 24.0; y = 6.0;
  printf("sm_pow(%i/1000, %i/1000)=%i", (int)(1000.0*x), (int)(y), (int)sm_pow(x, y));
  printf("sm_powint(%i/1000, %i/1000)=%i", (int)(1000.0*x), (int)(y), (int)sm_powint(x, y));
  x = 2.0; y = 0.5;
  printf("sm_pow(%i/1000, %i/1000)=%i/1000", (int)(1000.0*x), (int)(1000.0*y), (int)(1000.0*sm_pow(x, y)));
}
//-------------------------------------------------------------------------------------------------------- 
void deinit() {
  // Destroy Window
  window_destroy(main_window);
}
//-------------------------------------------------------------------------------------------------------- 
int main() {
  init();
  app_event_loop();
  deinit();
}
