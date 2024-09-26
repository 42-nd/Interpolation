# Interpolation Report
Цель: Сформировать практические навыки интерполяции табличных функций и численного дифференцирования.
Вариант номер 6.
![image](https://github.com/user-attachments/assets/7314220e-6053-49f8-bc97-16555d8d2395)

Расчет программы:
![image](https://github.com/user-attachments/assets/29756289-3e8e-48a7-be00-56f939f13bb3)
![image](https://github.com/user-attachments/assets/0173e7c8-77c1-49cd-937f-c02778645909)

![image](https://github.com/user-attachments/assets/047f91df-6909-4cf0-b7cc-26b51ae20f86)

![image](https://github.com/user-attachments/assets/c6db05ee-64ca-44e9-95c7-5527dd0b76e9)

![image](https://github.com/user-attachments/assets/15735caf-e5be-4f0d-a95f-5cd2f86fc120)

![image](https://github.com/user-attachments/assets/c8ed1db9-b044-4a42-ae99-d0d21dd304ef)
![image](https://github.com/user-attachments/assets/cb9a155c-1669-4b16-b93c-3796e38e024d)

По расчетам можно заметить, что для первой производной погрешность аппроксимации выше, чем для производной порядка 0, но ниже чем для второй производной.

Также хочется отметить 4 пункт лабораторной работы

![image](https://github.com/user-attachments/assets/bd2c99c0-03d1-4921-8b27-f66fe0cec437)

Как видно из последнего скриншота, значения отличаются, это можно объяснить тем, что формулы основываются на конечных разностях, что является линейным приближением первой производной. Этот метод приближён и учитывает только два соседних узла для оценки производной. Такой подход может давать неточные результаты.

Вывод:
Сплайн обеспечивает достаточно точную аппроксимацию функции и её первой производной.
Погрешности возрастают для второй производной, что характерно для методов интерполяции.
