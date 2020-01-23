// #include <iostream>
//
// template <typename Derived>
// class Base{
// public:
//   void interface(){
//     static_cast<Derived*>(this)->implementation();
//   }
//   void implementation(){
//     std::cout << "Implementation Basedf" << std::endl;
//   }
// };
//
// class Derived1: public Base<Derived1>{
// public:
//   void implementation(){
//     std::cout << "Implementation Derived1" << std::endl;
//   }
// };
//
// class Derived2: public Base<Derived2>{
// public:
//   void implementation(){
//     std::cout << "Implementation Derived2" << std::endl;
//   }
// };
//
// class Derived3: public Base<Derived3>{};
//
// template <typename T>
// void execute(T& base){
//   base.interface();
// }
//
// // arma::mat FisherScoring::GLMm{
// //
// //   std::cout << std::endl;
// //
// //   Derived1 d1;
// //   execute(d1);
// //
// //   Derived2 d2;
// //   execute(d2);
// //
// //   Derived3 d3;
// //   execute(d3);
// //
// //   std::cout << std::endl;
// // }
//
// RCPP_MODULE(BASEMODULE){
//   Rcpp::class_<Base>("Base")
//   .constructor()
//   .method("interface", &Base::interface)
//   ;
// }
//
// RCPP_MODULE(export_module_heritage){
//   Rcpp::class_<Base>("Base")
//   .constructor()
//   ;
//   Rcpp::class_<Derived2>("Derived2")
//     .derives<Base>("Base")
//     .constructor()
//     .method( "Derived2", &Derived2::implementation )
//   ;
// }
