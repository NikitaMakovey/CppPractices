# CppPractices
 
# mat_vector
For installing coverage you need make next steps:

   cd mat_vector/
   brew install lcov
   g++ --coverage main.cpp Vector.cpp Matrix.cpp -o main -std=c++17
   ./main
   lcov -t "main" -o main.info -c -d .
   genhtml -o report main.info
   
Open report/index.html into your Browser and check coverage.
 
 
