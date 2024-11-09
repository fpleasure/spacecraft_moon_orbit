# О проекте
Моделируется компланарное выведение КА с поверхности Луны на орбиту ИСЛ с заданными
параметрами. Рассматривается схема выведения, содержащая пассивный участок.

## Начало работы

### Содержание репозитория

```
├── img                   			# Графики
├── results                		 	# Результаты
├── src                    		    # Файлы программы
│   ├──files.py               		# Вывод в файл
│   ├──spacecraft_algorithm.py		# Алгоритм задачи
│   ├──parametrs.py					# Параметры
│   ├──runge_kutta.py               # Методы РК
│   ├──tests.py						# Тесты
│   ├──main.py						# Главный файл
├── requirements.txt
└── README.md
```

### Используемые технологии

- [Python](https://www.python.org/)

### Установка

Пошаговая инструкция для установки и запуска:

1. Клонируйте репозиторий:

    ```bash
    git clone https://github.com/fpleasure/spacecraft_moon_orbit.git
    ```

2. Перейдите в директорию проекта:

    ```bash
    cd spacecraft_moon_orbit
    ```

3. Создайте и активируйте виртуальное окружение:

    ```bash
    python -m venv venv
    source venv/bin/activate  # для Windows используйте `venv\Scripts\activate`
    ```

4. Установите необходимые зависимости:

    ```bash
    pip install -r requirements.txt
    ```

### Использование
Для запуска алгоритма используется функция **spacecraft_algorithm** из файла *spacecraft_algorithm*.

![Получившаяся траектория](/img/for_git.png)

## Контакты
- [Никита Королев](https://t.me/niki_korolev)
