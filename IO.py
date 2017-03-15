from reportlab.pdfgen import canvas;
from reportlab.lib.pagesizes import letter;
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

def pdf_save(data,file_name):
	Story=[];
	doc = SimpleDocTemplate(file_name, pagesize=letter,
                        rightMargin=72, leftMargin=72,
                        topMargin=72, bottomMargin=18);
	content = [];
	styles = getSampleStyleSheet()
	styles.fontsize = 18;
	styles.fontName = 'Helvetica'
	ptext = data;
	paragraphStyle = Paragraph(ptext, styles["Normal"]);
	Story.append(paragraphStyle);
	doc.build(Story);
		
