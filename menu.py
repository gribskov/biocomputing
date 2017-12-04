class Menu(object):
    '''
    combines a menu definition and action dispatcher
    '''
    keyseparator = ')'
    optionseparator = '  '

    def __init__(self):
        '''
        initialize empty menu object
        usage
            menu = Menu()
        '''
        self.all = {}
        self.query = {}
        self.response = ''

        return None

    def add(self, id, fields):
        '''
        add a menu to the object
        :fields: dictionary with the keys and strings for each option
        :return: True
        '''
        newmenu = {}
        query = ''
        for k in fields:
            newmenu[k] = fields[k].replace(k, '{0}{1}'.format(k.upper(), Menu.keyseparator))
        self.all[id] = newmenu

        query = Menu.optionseparator.join(newmenu.values())
        self.query[id] = query
        # print('id:', id, 'query:', query)

        return True

    def ask(self, id, indent):
        '''
        display the menu and get the response
        :param id: name of the menu
        :param indent: number of spaces to indent the query string
        :return: character read from keyboard
        usage
            while True:
                response = menu.ask('top', 2)
        '''
        while True:
            response = input('\n{0}{1}: '.format(' ' * indent, self.query[id]))
            response = response.upper()[0]
            try:
                if response in self.all[id]:
                    break
                else
                    print('Menu:{0} unknown option{1}\n'.format(id, response))
                    continue
            except:
                print('Menu:{0} unknown option{1}\n'.format(id, response))
                continue

        self.response = response
        return response


if __name__ == '__main__':
    print('Menu::Testing ')

    menu = Menu()
    menu.add('top', {'A': 'Add', 'Q': 'Quit'})
    response = menu.ask('top', 2)
    print('response is', response)