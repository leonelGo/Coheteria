from PIL import Image, ImageDraw, ImageFont

# Carga la imagen
image = Image.open('Logo_Altum_Principia.jpg')  # Asegúrate de que la ruta a la imagen sea correcta

# Define un contexto de dibujo
draw = ImageDraw.Draw(image)

# Define la posición y tamaño de la caja de texto para sobrescribir el texto incorrecto
text_position = (370, 630)  # Ajusta estos valores según sea necesario
text_box_size = (250, 50)  # Ajusta estos valores según sea necesario
background_color = (23, 32, 42)  # Ajusta este valor según sea necesario

# Cubre el texto antiguo con un rectángulo que coincida con el color de fondo
rectangle_upper_left = text_position
rectangle_lower_right = (text_position[0] + text_box_size[0], text_position[1] + text_box_size[1])
draw.rectangle([rectangle_upper_left, rectangle_lower_right], fill=background_color)

# Define el texto corregido y la fuente
corrected_text = "PRINCIPIA"
font = ImageFont.load_default()  # O utiliza ImageFont.truetype para una fuente específica y tamaño

# Calcula la posición para centrar el texto en la caja original del texto

text_x = text_position[0] + (text_box_size[0] - 100) / 2
text_y = text_position[1] + (text_box_size[1] - 100) / 2

# Dibuja el texto corregido en la imagen
draw.text((text_x, text_y), corrected_text, fill='#FF5733', font=font)  # Ajusta el color según sea necesario

# Guarda la imagen editada
edited_image_path = 'edited_image_corrected.jpg'
image.save(edited_image_path)

# Imprime la ruta de la imagen editada
print(edited_image_path)
