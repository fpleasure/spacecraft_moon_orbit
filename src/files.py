def initialize_output_file(filename="results/output.txt"):
    """
    Инициализирует файл для записи значений вектора u и времени t, создавая заголовок.

    Параметры:
    ----------
    - filename (str): Имя файла, в который будут записываться значения. По умолчанию "output.txt".

    Результат:
    ----------
    Создает или очищает файл с заданным именем и записывает в него заголовок для значений t и u.
    """
    try:
        with open(filename, "w") as file:
            file.write("t v_x v_y x y m\n")
        print(f"[INFO] Файл {filename} успешно инициализирован.")
    except Exception as e:
        print(f"[ERROR] Не удалось инициализировать файл {filename}: {e}")


def log_iteration_values(t, u, filename="results/output.txt"):
    """
    Записывает значения времени t и вектора u в файл на каждой итерации.

    Параметры:
    ----------
    - t (float): Текущее время итерации.
    - u (list of float): Вектор значений состояния системы на текущей итерации.
    - filename (str): Имя файла, в который будут записываться значения. По умолчанию "output.txt".

    Результат:
    ----------
    Добавляет строку с текущими значениями t и u в файл.
    """
    try:
        with open(filename, "a") as file:
            # Преобразуем значения времени и вектора u в строку
            u_str = " ".join(f"{val:.7f}" for val in u)
            line = f"{t:.7f} {u_str}\n"
            file.write(line)
    except Exception as e:
        print(f"[ERROR] Не удалось записать значения t и u в файл {
              filename}: {e}")
